#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <kseq.h>
#include <zlib.h>
#include <omp.h>
#include <unistd.h>
#include <argp.h>
#include <sys/types.h>
#include <sys/stat.h>
KSEQ_INIT(gzFile, gzread);

#define MATCH 0  /* enumerated type symbol for match */
#define INSERT 1 /* enumerated type symbol for insert */
#define DELETE 2 /* enumerated type symbol for delete */
const unsigned long int MAX_SIZE_EDIT_DIST_MATRIX = 512;

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define MAX2(a, b) ((a) > (b) ? (a) : (b))
#define MIN2(a, b) ((a) < (b) ? (a) : (b))

typedef struct
{
    int cost;   /* cost of reaching this cell */
    int parent; /* parent cell */
} cell;

/*--------------- Functions prototypes------------------------*/

/* Function in which a DNA sequence has been processed, given a FASTA formatted file.
 * If more sequences are present in the input file the function will be recalled recursively until the end of it.*/
void sequences_reader(char *pattern, kseq_t *seq, double k_edit_dist_perc, int num_threads);

/* Functions that return the length for DNA subsequence, given the pattern and whole sequence itself.*/
unsigned int subdivide_sequence(char *pattern, kseq_t *sequence);

/* Function that compute the edit distance matrix. This needs as inputs both previous and current matrix
 * in order to correctly initialize the current matrix. Boolean show_mat it has been used for debugging purpose*/
void string_compare_parallel(char *p, char *t, int nr, int nc, int ncol_prev, cell *mat, cell *mat_p ,bool show_mat, int num_threads);

/* Function that it is responsible of finding the j-th indexes in the current matrix.
 * Return the total number of matches for the current matrix.*/
int *k_dist_submatch_indexes(char *s, char *t,int *num_ind_k, int k_edit_dist, int nc, cell *mp, cell *m);

/* Function that initialize the first column of edit distance matrix in case of first sub matrix.*/
void column_init_submatch(int nr, int nc, cell *matrix);

/* Function that initialize the first column of edit distance matrix in case of sub matrix different from the first one.*/
void column_init_submatch_next_win(int nr, int nc, cell *matrix_curr, int ncol_prev, cell *matrix_prev, char *pattern, char *t);

/* Function that initialize the first row of edit distance matrix*/
void row_init_submatch(int nr, int nc, cell *matrix);

/*Function that will save the edit path along with related substrings, saving these results in substrings matrix passed in input.*/
void susbtrings_matched(char *p, char *t, char *text_prev, int *indexes, int n_submatch, int nr, int nc, cell *mat,
                        int max_rows_sub, int max_col_sub, cell *matp, int ncol_prev, char *subtrings);

/* Save on disk edit paths and substrings found
 * In this file also the pattern, the k threshold and length of current subsequence has been passed in order to save the results in appropriate files
 * for different cases.*/
void save_subtrings_matched(char *p, int edit_dist, unsigned long int length_seq, int nr, int nc, char *substrings, cell *matp);

/*Function that will be recalled recursively and that will find the edit path, given as input an index j and both current and previous sub matrix*/
void reconstruct_path(char *s, char *t, char *t_prev,int i, int j, char *edit_path, int nr, int nc, cell *matp, int ncol_prev, cell *mat);

/* match_out, delete_out, insert_out are all functions used for building the string of edit path that will be saved as output*/
void match_out(char *s, char *t, char *t_prev, int i, int j, cell *mat_prev, char *edit_path);
void delete_out(char *s, int i, char *edit_path);
void insert_out(char *t, int j, char *edit_path);

/*Function that given and edit path will find the correspondent substring*/
void substring_from_text(char *sub, char *text, char *text_prev, cell *matp ,char *edit_path, int last_char_ind);

/*-----------------------------------------------------------*/

int main(int argc, char *argv[])
{
    //Input data
	gzFile fp_pattern, fp;
    kseq_t *seq, *seq_pattern;
    int l, l_pattern;
    char* pattern;
    char* file_pattern_name = argv[1];
    char* file_sequences_name = argv[2];
    int k_edit_dist_perc = atoi(argv[3]);
    int num_threads = atoi(argv[4]);

    // Open pattern file
    fp_pattern = gzopen(file_pattern_name,"r");
    seq_pattern = kseq_init(fp_pattern);

    struct stat st = {0};

    // Create outputs directory if not exist
    if (stat("./outputs", &st) == -1) {
    	mkdir("./outputs", 0700);
    }

    // This while loop will compare each pattern found in pattern file with each DNA sequence found in FASTA input file
    while((l_pattern = kseq_read(seq_pattern)) >= 0){
    	pattern = calloc((strlen(seq_pattern->seq.s)+1), sizeof(char));
    	pattern = strcpy(pattern, seq_pattern->seq.s);

    	char dir_name_pattern_results[100]; // 100 Length filename in which save results, including also the character terminator
    	sprintf(dir_name_pattern_results, "./outputs/results_pattern[%li]", strlen(pattern));

    	// Create results directory for the current pattern length
    	if (stat(dir_name_pattern_results, &st) == -1) {
    	    mkdir(dir_name_pattern_results, 0700);
    	}

    	fp = gzopen(file_sequences_name, "r");
    	seq= kseq_init(fp);


    	if ((l = kseq_read(seq)) >= 0){
    		sequences_reader(pattern, seq, k_edit_dist_perc, num_threads);
    	}

    	free(pattern);
    	kseq_destroy(seq);
    	gzclose(fp);

    }

    kseq_destroy(seq_pattern);
    gzclose(fp_pattern);


    return 0;
}

/*----------------------------------Functions definition --------------------------------------*/

void sequences_reader(char *pattern, kseq_t *seq, double k_edit_dist_perc, int num_threads){
    char *name_seq;
    int l;
    int num_win = 0;
    int k_edit_dist = (k_edit_dist_perc/100)*strlen(pattern);
    cell *Mcurr, *Mprev;
    char  *text_curr, *text_prev;
    int ncol_curr, ncol_prev;
    int nrows = strlen(pattern) + 1;
    unsigned long int length_sequence = strlen(seq->seq.s);
    bool show = false;

    /* Create directory for result given a certain sequence and pattern length */
    struct stat st = {0};
    char dir_name_results_text[100];
	sprintf(dir_name_results_text, "./outputs/results_pattern[%li]/text[%li]",strlen(pattern),length_sequence);
	if (stat(dir_name_results_text, &st) == -1) {
		mkdir(dir_name_results_text, 0700);
	}
	/*--------------------------------------------------------------*/
	// This do while loop iterate until we are in the same DNA sequence
	double initialTime_searching_sequence = omp_get_wtime();
    do{
    	name_seq = calloc((strlen(seq->name.s)+1), sizeof(char));
    	name_seq = strcpy(name_seq, seq->name.s);
        /* Dividing sequence in subsequences in order to use a bounded size of memory */
        int num_char_subsequence = subdivide_sequence(pattern, seq);
        for(int index_char=0; index_char < length_sequence; index_char += num_char_subsequence){

        	/* Get subsequence from the complete sequence string */
        	if(index_char + num_char_subsequence >= length_sequence){
        		num_char_subsequence = length_sequence - index_char;
        	}
        	char *subsequence = calloc(num_char_subsequence+1,sizeof(char));
        	memcpy(subsequence, &seq->seq.s[index_char], num_char_subsequence);
        	subsequence[num_char_subsequence] = '\0';
        	/*-------------------------------------------------*/
        	/* Use the subsequence to create the relative edit distance matrix with the pattern
        	 * num_win variable will count the number of sub matrix that have to be computed, and has been used only for making a different initialization for the current and previous matrix */
        	num_win++;
			if (num_win == 1) {
				text_curr = subsequence;
				ncol_curr = strlen(text_curr) + 1;
				Mcurr = calloc(nrows * ncol_curr, sizeof(cell));
				Mprev = NULL;
				text_prev = NULL;
			} else {
				if (num_win > 2) {
					free(text_prev);
					text_prev = NULL;
				}
				text_prev = text_curr;
				ncol_prev = ncol_curr;
				if (num_win > 2) {
					free(Mprev);
					Mprev = NULL;
				}
				Mprev = Mcurr;
				Mcurr = NULL;
				text_curr = NULL;
				text_curr = subsequence;
				ncol_curr = strlen(text_curr);
				Mcurr = calloc(nrows * ncol_curr, sizeof(cell));
			}

			// instantiation of dynamic programming table for edit distances
			string_compare_parallel(pattern, text_curr, nrows, ncol_curr, ncol_prev, Mcurr, Mprev, show, num_threads);

			int num_submatch;
			int *submatch_indexes;

			// Indexes searching
			submatch_indexes = k_dist_submatch_indexes(pattern, text_curr, &num_submatch, k_edit_dist, ncol_curr, Mprev, Mcurr);

			double initialTime_substrings_searching = omp_get_wtime();

            #pragma omp parallel firstprivate(k_edit_dist, length_sequence, ncol_curr, ncol_prev, nrows) num_threads(num_threads)
			{

				int *sub_indexes_k_thread;
				int num_ind_k_thread = 0;

				// Here the total set of j-th indices has been spread equally among the threads
                #pragma omp for schedule(static)
				for (int k = 0; k < num_submatch; k++) {
					if (num_ind_k_thread == 0) {
						num_ind_k_thread = num_ind_k_thread + 1;
						sub_indexes_k_thread = malloc(sizeof(int));
						sub_indexes_k_thread[num_ind_k_thread - 1] = submatch_indexes[k];
					} else {
						num_ind_k_thread = num_ind_k_thread + 1;
						sub_indexes_k_thread = realloc(sub_indexes_k_thread, num_ind_k_thread * sizeof(int));
						sub_indexes_k_thread[num_ind_k_thread - 1] = submatch_indexes[k];
					}
				}

				int max_lenght_substring_matched_t = strlen(pattern) + k_edit_dist + 1;
				int max_rows_thread = 2 * num_ind_k_thread; // Because  each pair of edit path and substring will be saved
				char *substring_matched_thread;
				substring_matched_thread = calloc(max_rows_thread * max_lenght_substring_matched_t, sizeof(char));

				// Regarding the subset of the total indices, each thread will start in searching the relative edit path along with substring
				susbtrings_matched(pattern, text_curr, text_prev,
						sub_indexes_k_thread, num_ind_k_thread, nrows,
						ncol_curr, Mcurr, max_rows_thread,
						max_lenght_substring_matched_t, Mprev, ncol_prev,
						substring_matched_thread);

				//Each thread will save the results in different files
                save_subtrings_matched(pattern, k_edit_dist, length_sequence, max_rows_thread, max_lenght_substring_matched_t, substring_matched_thread, Mprev);
                free(substring_matched_thread);
			}

			double finalTime_substrings_searching = omp_get_wtime();
	    	printf("\nElapsed time total for substrings searching: %f [s]\n\n",finalTime_substrings_searching - initialTime_substrings_searching);

        }

        l = kseq_read(seq);

    } while ((l >= 0) & (seq->name.s == name_seq));

    double finalTime_searching_sequence = omp_get_wtime();
	printf("\n\nELAPSED TOTAL TIME SEQUENCE: %f [s]\nPATTERN LENGTH [%li]\nSEQ LENGTH[%li]\nK_EDIT_DIST[%i]\nTHREAD NUMBERS:[%i]\n\n",
			finalTime_searching_sequence - initialTime_searching_sequence,
			strlen(pattern), length_sequence, k_edit_dist, num_threads);

	FILE *fptr_results;
	fptr_results = fopen("./results_time.txt", "a");
	fprintf(fptr_results,"\n\nELAPSED TOTAL TIME SEQUENCE: %f [s]\nPATTERN LENGTH [%li]\nSEQ LENGTH[%li]\nK_EDIT_DIST[%i]\nTHREAD NUMBERS:[%i]\n\n",
			finalTime_searching_sequence - initialTime_searching_sequence, strlen(pattern), length_sequence, k_edit_dist, num_threads);
	fclose(fptr_results);

    if (l == -1){
    	free(Mcurr);
    	free(Mprev);
    	free(text_curr);
    	free(text_prev);
        return;
    }else if (seq->name.s != name_seq){
    	free(Mcurr);
    	free(Mprev);
    	free(text_curr);
    	free(text_prev);
        sequences_reader(pattern, seq, k_edit_dist_perc, num_threads);
    }

    return;
}


void string_compare_parallel(char *p, char *t, int nr, int nc, int ncol_prev, cell *mat, cell *mat_p ,bool show_mat, int num_threads)
{

	int opt[3]; /* cost of the three options */
	int plen = strlen(p);
	int length_curr_diag;
	int num_diagonals = nc + nr - 1;
	int first_col, first_row;

	row_init_submatch(nr, nc, mat);
	// This check will change the initialization of the first column in relation if we are computing the first sub matrix or not
	if (mat_p == NULL) {
		column_init_submatch(nr, nc, mat);
	} else {
		column_init_submatch_next_win(nr, nc, mat, ncol_prev, mat_p, p, t);
	}

	double initialTime_edit_distance_matrix = omp_get_wtime();
	for (int c = 2; c < num_diagonals; c++) {

		length_curr_diag = MIN2(plen + 1, c + 1);
		first_col = MIN2(c - 1, nc - 1);
		first_row = MAX2(1, c - nc + 1);

		int n_iterations = length_curr_diag - first_row;

		// This control has to be made in order to not compute element on diagonal that belongs to the first column, that is already initialized
		if (first_col - n_iterations < 0) {
			n_iterations--;
		}

		// This offset it's necessary when a sub matrix different from the first one has to be computed
		// so the index for comparing pattern and subsequence of DNA are correctly indexed.
		int offset_text_index;
		(mat_p == NULL) ? (offset_text_index = -1) : (offset_text_index = 0);

		// This scaling has to be performed when the number of threads requested are equal to the maximum available
		int num_max_threads = omp_get_max_threads();
		if(num_threads == num_max_threads){
			num_threads = num_threads -2;
		}else if(num_threads > num_max_threads){
			num_threads = num_max_threads -2;
		}


		#pragma omp parallel for num_threads(num_threads) firstprivate(nc, first_col, first_row, n_iterations, offset_text_index) private(opt)
		for (int col = first_col; col > first_col - n_iterations; col--) {
			int row = first_row + (first_col - col);
			int mat_index = nc*row + col;
			int match_cost = (p[row-1] == t[col + offset_text_index]) ? 0 : 1;

			opt[MATCH] = mat[mat_index - nc - 1].cost + match_cost;
			opt[INSERT] = mat[mat_index - 1].cost + 1;
			opt[DELETE] = mat[mat_index - nc].cost + 1;
			mat[mat_index].cost = opt[MATCH];
			mat[mat_index].parent = MATCH;
			for (int k = INSERT; k <= DELETE; k++){
				if (opt[k] < mat[mat_index].cost) {
					mat[mat_index].cost = opt[k];
					mat[mat_index].parent = k;
				}
			}
		}
	}

	double finalTime_edit_distance_matrix = omp_get_wtime();
	printf("\nElapsed time for edit_distance matrix: %f [s]\n\n",finalTime_edit_distance_matrix-initialTime_edit_distance_matrix);

	// Only used for debug purposes
	if(show_mat == true) print_matrix(p, t, nr, nc, mat);
}


int *k_dist_submatch_indexes(char *s, char *t,int *num_ind_k, int k_edit_dist, int nc, cell *mp, cell *m)
{
    int i, offset;
    // Offset serves because indices are saved considered as indices of the current DNA subsequence
    // So if first sub matrix has to be computed then it's necessary to subtract 1 otherwise it's
    // possible to take index with its value
    if (mp == NULL){
        i = 1;
        offset = -1;
    }else{
        i = 0;
        offset = 0;
    }

    int last_row = strlen(s);
    int *sub_indexes_k;
    *num_ind_k = 0;

	for (int k = i; k < nc; k++) {
		int mat_index = nc * last_row + k;
		if (m[mat_index].cost <= k_edit_dist) {
			if (*num_ind_k == 0) {
				*num_ind_k = *num_ind_k + 1;
				sub_indexes_k = malloc(sizeof(int));
				sub_indexes_k[*num_ind_k - 1] = k + offset;
			} else {
				*num_ind_k = *num_ind_k + 1;
				sub_indexes_k = realloc(sub_indexes_k,
						*num_ind_k * sizeof(int));
				sub_indexes_k[*num_ind_k - 1] = k + offset;
			}
		}
	}



    return sub_indexes_k;
}


void susbtrings_matched(char *p, char *t, char *text_prev, int *indexes, int n_submatch, int nr, int nc, cell *mat,
                        int max_rows_sub, int max_col_sub, cell *matp, int ncol_prev, char *subtrings)
{

    int last_row = strlen(p);
    char edit_path[max_col_sub];
    char substring[max_col_sub];
    memset(edit_path, 0, max_col_sub);
    memset(substring, 0, max_col_sub);

    int offset;
    // Depending on we are referring to the first or not sub matrix indices has to be shifted in order to correctly reach characters in DNA subsequence
    (matp == NULL) ? (offset = 1)  : (offset = 0);

    for (int j = 0; j < n_submatch; j++)
    {
        int index_last_char_sub = indexes[j] + offset;
        reconstruct_path(p, t, text_prev, last_row, index_last_char_sub, edit_path, nr, nc, matp, ncol_prev, mat);
        substring_from_text(substring, t, text_prev, matp, edit_path, index_last_char_sub - offset);

        int index_sub_mat;
        for (int c = 0; c < max_col_sub; c++){
        	// This serves in order to correctly index the substrings matrix: given an edit path and a substring
        	// these strings are saved one row at time, so for a new pair it's necessary to move forward of two rows.
            index_sub_mat = max_col_sub*(j*2) + c;
             //Row for edith path
            subtrings[index_sub_mat] = edit_path[c];
             //Row for text substring
            index_sub_mat = index_sub_mat + max_col_sub;
            subtrings[index_sub_mat] = substring[c];
        }

        memset(edit_path, 0, max_col_sub);
        memset(substring, 0, max_col_sub);

    }
}


void save_subtrings_matched(char *p, int edit_dist, unsigned long int length_seq ,int nr, int nc, char *substrings, cell *matp){

	int p_length = strlen(p);
	int thread_num = omp_get_thread_num();
	char file_name_saved_results[150]; // 150 Length filename in which save results, including also the character terminator
	sprintf(file_name_saved_results, "./outputs/results_pattern[%li]/text[%li]/susbtrings_matched_thread[%i]_kdist[%i].txt.gz", p_length,length_seq,thread_num,edit_dist);

	gzFile fptr;
    if (matp == NULL){
        fptr = gzopen(file_name_saved_results,"w");
        gzprintf(fptr, "\n---LEGEND---\n");
        gzprintf(fptr, "\n[I]:Insert operation\n[D]:Delete operation\n[S]:Sobstitution operation\n[M]:Match");
        gzprintf(fptr, "\nAll operations cost 1 except for [M]atch that does not involve any operation.");
        gzprintf(fptr, "\n\n------------\n");
        gzprintf(fptr, "\nPATTERN: %s\n", p);
        gzprintf(fptr, "\nEdist distance: %d\n\n", edit_dist);
    }else{
        fptr = gzopen(file_name_saved_results,"a");
    }

    for (int r = 0; r < nr; r++){

    	// This is uses to print the correct number for the related match
    	// Each pair of edit path and substring will correspond to one match, and so the number of match has been set
        if (r == 0){
        	gzprintf(fptr,"[%d]\n",r+1);
        }else if ( r%2 == 0){
            int num_occurence = (r/2)+1;
            gzprintf(fptr,"\n[%d]\n", num_occurence);
        }

        int c=0;
        int index_char=nc*r+c;
        do{
        	gzprintf(fptr,"%c", substrings[index_char]);
        	index_char++;
        }while(substrings[index_char] != '\0');

        gzprintf(fptr,"\n\n");
    }
    gzclose(fptr);
}

unsigned int subdivide_sequence(char *pattern, kseq_t *sequence){

	unsigned long int max_size_matrix = MAX_SIZE_EDIT_DIST_MATRIX*1048576; // Dimension expressed in Byte
	unsigned long int length_sequence = strlen(sequence->seq.s)*2; // Multiplied by two in order to check if the length of current sequence could be handled entirely
	unsigned long int length_pattern = strlen(pattern);
	unsigned long int length_subsequence =  length_sequence, size_matrix_subseqeuence;
	unsigned long int size_cell = sizeof(cell);
	do{
		length_subsequence = length_subsequence/2;
		size_matrix_subseqeuence = (length_subsequence * length_pattern * size_cell);
	}while(size_matrix_subseqeuence > max_size_matrix);

	return length_subsequence;
}

/*-------------------------Functions useful for edit distance matrix initialization-------------------------*/
void row_init_submatch(int nr, int nc, cell *matrix)
{
    for (int h = 0; h < nc; h++)
    {
        matrix[h].cost = 0;
        matrix[h].parent = -1; /* note change */
    }
}

void column_init_submatch(int nr, int nc, cell *matrix)
{
    for (int h = 0; h < nr; h++)
    {
        matrix[nc*h].cost = h;
        if (h > 0)
            matrix[nc*h].parent = DELETE;
        else
            matrix[nc*h].parent = -1;
    }
}

void column_init_submatch_next_win(int nr, int nc, cell *matrix_curr, int ncol_prev, cell *matrix_prev, char *pattern, char *t){
    int opt[3];

    // This function will initialize the first column of a new sub matrix in case which this one it is not the first one.
    // So the first column has to be initialized taking into account the values of last column of the previous sub matrix
    // Evaluating edit costs and type of edit operations accordingly to edit distance algorithm

    for (int h = 1; h < nr; h++)
    {

        int mat_curr_index = nc*h;
        int mat_prev_index_diag = ncol_prev*(h-1) + (ncol_prev-1);
        int mat_prev_index_row_move = ncol_prev*h + (ncol_prev-1);
        int mat_curr_index_col_move = nc*(h-1);
        int match_cost = (pattern[h - 1] == t[0]) ? 0 : 1;
        opt[MATCH] = matrix_prev[mat_prev_index_diag].cost + match_cost;
        opt[INSERT] = matrix_prev[mat_prev_index_row_move].cost + 1;
        opt[DELETE] = matrix_curr[mat_curr_index_col_move].cost + 1;
        matrix_curr[mat_curr_index].cost = opt[MATCH];
        matrix_curr[mat_curr_index].parent = MATCH;
        for (int k = INSERT; k <= DELETE; k++)
            if (opt[k] < matrix_curr[mat_curr_index].cost)
            {
                matrix_curr[mat_curr_index].cost = opt[k];
                matrix_curr[mat_curr_index].parent = k;
            }
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------Functions useful for reconstructing edit operations and for finding substrings ----------------------*/
void insert_out(char *t, int j, char *edit_path)
{
    char ins = 'I';
    strncat(edit_path, &ins, 1);
}

void delete_out(char *s, int i, char *edit_path)
{
    char del = 'D';
    strncat(edit_path, &del, 1);
}

void match_out(char *s, char *t, char *t_prev, int i, int j, cell *mat_prev, char *edit_path)
{
    int index_text_char;
    char match = 'M';
    char mismatch = 'S';

    if (j < 0){
        t = t_prev;
        index_text_char = strlen(t_prev) + j;
    }else{
        (mat_prev == NULL) ? (index_text_char = j-1) : (index_text_char = j);
    }

    if (s[i - 1] == t[index_text_char])
        strncat(edit_path, &match, 1);
    else
        strncat(edit_path, &mismatch, 1);
}

void reconstruct_path(char *s, char *t, char *t_prev,int i, int j, char *edit_path, int nr, int nc, cell *matp, int ncol_prev, cell *mat)
{
    int mat_index;
    /*This is necessary for dealing with possible substrings matched at the boundary between two sub matrices
    * So it's necessary to change the reference passing to previous matrix and accordingly change how index the matrix.*/
    if (j < 0){
        mat = matp;
        mat_index = ncol_prev*i + (ncol_prev + j);
    }else{
        mat_index = nc*i +j;
    }


    if (mat[mat_index].parent == -1)
        return;
    if (mat[mat_index].parent == MATCH)
    {
        reconstruct_path(s, t, t_prev, i - 1, j - 1, edit_path, nr, nc, matp, ncol_prev, mat);
        match_out(s, t, t_prev, i, j, matp, edit_path);
        return;
    }
    if (mat[mat_index].parent == INSERT)
    {
        reconstruct_path(s, t, t_prev, i, j - 1, edit_path, nr, nc, matp, ncol_prev, mat);
        insert_out(t, j, edit_path);
        return;
    }
    if (mat[mat_index].parent == DELETE)
    {
        reconstruct_path(s, t, t_prev, i - 1, j, edit_path, nr, nc, matp, ncol_prev, mat);
        delete_out(s, i, edit_path);
        return;
    }
}



void substring_from_text(char *sub, char *text, char *text_prev, cell *matp ,char *edit_path, int last_char_ind){

	// Given the j index of the match will move backwards on a number of positions equal to length of edit path
	int first_char = last_char_ind - (strlen(edit_path) - 1);
    if (first_char < 0){
         if (matp == NULL){
        	 // If going backwards led to negative index for first char has been found but we are in the first sub matrix
        	 // then the first char is set to zero and the related substring has been retrieve
             first_char = 0;
             for (int i = first_char; i <= last_char_ind; i++) sub[i-first_char] = text[i];
         }else{
        	 // If going backwards negative index for first char has been found but we are not in the first sub matrix
        	 // then the substring has been retrieved as a part form the current DNA subsequence and a part form the previous subsequence
             int first_char_curr_text = 0;
             int last_char_prev_text = strlen(text_prev) - 1;
             int first_char_prev_text = last_char_prev_text + (first_char + 1);
             int i = 0;
             for (int s = first_char_prev_text; s <= last_char_prev_text; s++){
                 sub[i] = text_prev[s];
                 i++;
             }
             for (int d = first_char_curr_text; d <= last_char_ind; d++){
                 sub[i] = text[d];
                 i++;
             }
         }
    }else{
        for (int i = first_char; i <= last_char_ind; i++) sub[i-first_char] = text[i];
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/


/*Used only for debugging purposes*/
void print_matrix(char *s, char *t, int nr, int nc, cell *mat)
{
    printf("\n      ");
    for (int i = 0; i < strlen(t); i++)
        printf("[%c]", t[i]);
    printf("\n");

    for (int row = 0; row < nr; row++)
    {
        (row == 0) ? printf("   ") : printf("[%c]", s[row - 1]);
        for (int column = 0; column < nc; column++)
        {
            int mat_index = nc*row + column;
            printf(" %d ", mat[mat_index].parent);
        }
        printf("\n");
    }
}






