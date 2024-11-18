//# Copyright Cristina Iordache 314CAa 2023-2024
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//tartarea overflowlui si numerele negative
int mod(int x)
{
	x = x % 10007;
	if (x < 0)
		x = x + 10007;
	return x;
}

//dubleaza capacitatea
int needs_extend(int *capacity, int *length)
{
	if (*capacity == 0 || *length == *capacity) {
		if (*capacity == 0)
			*capacity = 1;
		else
			*capacity = 2 * (*capacity);
		return 1;
	}
	return 0;
}

//functie pentru alocare memorie matrice
int **alloc_matrix(int m, int n)
{
	int **mat = (int **)malloc(m * sizeof(int *));
	if (!mat)
		exit(EXIT_FAILURE);
	for (int i = 0; i < m; i++) {
		mat[i] = (int *)calloc(n, sizeof(int));
		if (!mat[i]) {
			for (int j = 0; j < i; j++)
				free(mat[j]);
			free(mat);
			exit(EXIT_FAILURE);
		}
	}
	return mat;
}

//functie eliberare matrice
void free_m(int m, int **matrix)
{
	for (int i = 0; i < m; i++)
		free(matrix[i]);
	free(matrix);
}

//functie inmultire matrici
int **multiply(int **a, int **b, int m, int n)
{
	int **rez = alloc_matrix(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				rez[i][j] = mod(rez[i][j] + mod(a[i][k] * b[k][j]));
	return rez;
}

//micsoreaza capacitatea
int needs_shrenk(int *capacity, int *length)
{
	if (*length < *capacity / 2) {
		if (*capacity == 0)
			*capacity = 1;
		else
			*capacity = (*capacity) / 2;
		return 1;
	}
	return 0;
}

//mareste dimensiunea array-ului
int ***extend_matrices(int ***matrices, int *capacity)
{
	matrices = (int ***)realloc(matrices, sizeof(int **) * (*capacity));
	return matrices;
}

//mareste matricea de dimensiuni
int **extend_size_matrix(int **sizematrix, int *capacity)
{
	sizematrix = realloc(sizematrix, sizeof(int *) * (*capacity));
	return sizematrix;
}

//micsoreaza dimensiunea array-ului
int ***shrenk_matrix(int ***matrices, int *capacity)
{
	matrices = (int ***)realloc(matrices, sizeof(int **) * (*capacity));
	return matrices;
}

//micsoreaza matricea de dimensiuni
int **shrenk_size_matrix(int **sizematrix, int *capacity)
{
	sizematrix = realloc(sizematrix, sizeof(int *) * (*capacity));
	return sizematrix;
}

//se introde in array o noua matrice
void read_new_matrix(int ***matrices, int **sizematrix, int *length)
{
	int m, n, x;
	scanf("%d%d", &m, &n);

	sizematrix[*length] = (int *)malloc(sizeof(int) * 2);
	sizematrix[*length][0] = m;
	sizematrix[*length][1] = n;

	matrices[*length] = (int **)malloc(sizeof(int *) * m);
	for (int i = 0; i < m; i++)
		matrices[*length][i] = (int *)calloc(n, sizeof(int));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) {
			scanf("%d", &x);
			matrices[*length][i][j] = x;
		}
	*length += 1;
}

//redimensionarea unei matrici
void resize(int ****matrices, int ***sizematrix, int *length)
{
	int index;
	int num_l, num_c;
	int pos_l = 0, pos_c = 0;
	int *lines;
	int *columns;
	int **tmp;
	int ***copy = *matrices;
	int **size_copy = *sizematrix;
	scanf("%d", &index);
	scanf("%d", &num_l);

	lines = (int *)malloc(num_l * sizeof(int));
	for (int i = 0; i < num_l; i++)
		scanf("%d", &lines[i]);
	scanf("%d", &num_c);
	columns = (int *)malloc(num_c * sizeof(int));
	for (int i = 0; i < num_c; i++)
		scanf("%d", &columns[i]);

	if (index >= *length || index < 0) {
		printf("No matrix with the given index\n");
		free(lines);
		free(columns);
		return;
	}
	tmp = alloc_matrix(num_l, num_c);
	for (int i = 0; i < num_l; i++) {
		pos_l = lines[i];
		for (int j = 0; j < num_c; j++) {
			pos_c = columns[j];
			tmp[i][j] = copy[index][pos_l][pos_c];
		}
	}
	free_m(size_copy[index][0], copy[index]);
	copy[index] = tmp;

	size_copy[index][0] = num_l;
	size_copy[index][1] = num_c;

	free(lines);
	free(columns);
}

void dimension(int **sizematrix, int *length)
{
	int index;
	scanf("%d", &index);
	if (index >= *length || index < 0) {
		printf("No matrix with the given index\n");
		return;
	}
	printf("%d %d\n", sizematrix[index][0], sizematrix[index][1]);
}

void print_matrix(int ***matrices, int **sizematrix, int *length)
{
	int index;
	scanf("%d", &index);
	if (index >= *length || index < 0) {
		printf("No matrix with the given index\n");
		return;
	}
	int m = sizematrix[index][0];
	int n = sizematrix[index][1];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%d ", matrices[index][i][j]);
		printf("\n");
	}
}

//inmultirea a doua matrici de index diferit
void mat_multiplication(int ****matrices, int ***sizematrix, int *length)
{
	int ***copy = *matrices;
	int **size_copy = *sizematrix;
	int index1, index2;
	int **tmp;
	scanf("%d%d", &index1, &index2);
	if (index1 >= *length || index2 >= *length || index1 < 0 || index2 < 0) {
		printf("No matrix with the given index\n");
		return;
	}
	int m1 = size_copy[index1][0];
	int n1 = size_copy[index1][1];
	int m2 = size_copy[index2][0];
	int n2 = size_copy[index2][1];
	if (n1 != m2) {
		printf("Cannot perform matrix multiplication\n");
		return;
	}
	tmp = alloc_matrix(m1, n2);
	for (int i = 0 ; i < m1; i++)
		for (int j = 0; j < n2; j++)
			for (int k = 0; k < m2; k++)
				tmp[i][j] = mod(tmp[i][j]
							+ mod(copy[index1][i][k] * copy[index2][k][j]));
	size_copy[*length] = (int *)malloc(sizeof(int) * 2);
	size_copy[*length][0] = m1;
	size_copy[*length][1] = n2;
	copy[*length] = tmp;
	*length += 1;
}

//transpusa unei matrici
void transpose(int ****matrices, int ***sizematrix, int *length)
{
	int ***copy = *matrices;
	int **size_copy = *sizematrix;
	int index;
	scanf("%d", &index);
	if (index >= *length || index < 0) {
		printf("No matrix with the given index\n");
		return;
	}
	int **tmp;
	int m = size_copy[index][0];
	int n = size_copy[index][1];
	tmp = alloc_matrix(n, m);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			tmp[i][j] = copy[index][j][i];

	free_m(m, copy[index]);
	copy[index] = tmp;

	int aux;
	aux = size_copy[index][0];
	size_copy[index][0] = size_copy[index][1];
	size_copy[index][1] = aux;
}

//suma elementelor unei matrici
int sum(int **mat, int m, int n)
{
	int s = 0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			s = mod(s + mat[i][j]);
	return s;
}

//sortarea array-ului de matrici
void matrix_sort(int ****matrices, int ***sizematrix, int *length)
{
	int ***copy = *matrices;
	int **size_copy = *sizematrix;
	int **aux;
	int aux_dim_r;
	int aux_dim_c;
	for (int i = 0; i < *length - 1; i++) {
		for (int j = i + 1; j < *length; j++) {
			int x = sum(copy[i], size_copy[i][0], size_copy[i][1]);
			int y = sum(copy[j], size_copy[j][0], size_copy[j][1]);
			if (x > y) {
				aux = copy[i];
				copy[i] = copy[j];
				copy[j] = aux;

				aux_dim_r = size_copy[i][0];
				size_copy[i][0] = size_copy[j][0];
				size_copy[j][0] = aux_dim_r;

				aux_dim_c = size_copy[i][1];
				size_copy[i][1] = size_copy[j][1];
				size_copy[j][1] = aux_dim_c;
			}
		}
	}
}

//ridicarea unei matrici la putere in timp logaritmic
void logarithmic_pow(int ****matrices, int ***sizematrix, int *length)
{
	int p;
	int index;
	int ***copy = *matrices;
	int **size_copy = *sizematrix;
	scanf("%d%d", &index, &p);
	if (index >= *length || index < 0) {
		printf("No matrix with the given index\n");
		return;
	}
	if (p < 0) {
		printf("Power should be positive\n");
		return;
	}
	int m = size_copy[index][0];
	int n = size_copy[index][1];
	if (m != n) {
		printf("Cannot perform matrix multiplication\n");
		return;
	}
	int **tmp;
	tmp = alloc_matrix(m, n);
	int **ans;
	ans = alloc_matrix(m, n);
	for (int i = 0; i < m ; i++) {
		for (int j = 0 ; j < n; j++) {
			tmp[i][j] = copy[index][i][j];
			ans[i][j] = (i == j) ? 1 : 0;
		}
	}
	int **newans;
	for (int bit = 0; bit < 31; bit++) {
		if (p & (1 << bit)) {
			newans = multiply(ans, tmp, m, n);
			free_m(m, ans);
			ans = newans;
		}
		newans = multiply(tmp, tmp, m, n);
		free_m(m, tmp);
		tmp = newans;
	}
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			copy[index][i][j] = ans[i][j];
	free_m(m, ans);
	free_m(m, tmp);
}

//elibearrea unei matrici de la un anumit index
void free_matrix(int ****matrices, int *length, int ***sizematrix)
{
	int index;
	int ***copy = *matrices;
	int **size_copy = *sizematrix;
	scanf("%d", &index);
	if (index < 0 || index >= *length) {
		printf("No matrix with the given index\n");
		return;
	}
	for (int i = 0; i < size_copy[index][0]; i++)
		free(copy[index][i]);
	free(copy[index]);
	free(size_copy[index]);

	int l = *length - 1;
	for (int i = index; i < l; i++) {
		copy[i] = copy[i + 1];
		size_copy[i] = size_copy[i + 1];
	}
	*length -= 1;
}

void plus(int **a, int **b, int **c, int m)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			c[i][j] = mod(a[i][j] + b[i][j]);
}

void minus(int **a, int **b, int **c, int m)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			c[i][j] = mod(a[i][j] - b[i][j]);
}

//inmultirea matricilor cu algoritmul Strassen
void strassen_multiply(int **a, int **b, int **c, int m)
{	int **a11, **a12, **a21, **a22;
	int **b11, **b12, **b21, **b22;
	int **c11, **c12, **c21, **c22;
	int **m1, **m2, **m3, **m4, **m5, **m6, **m7;
	int **aresult, **bresult;
	if (m == 1) {
		c[0][0] = a[0][0] * b[0][0];
	} else {
		int n = m / 2;
		a11 = alloc_matrix(n, n); a12 = alloc_matrix(n, n);
		a21 = alloc_matrix(n, n); a22 = alloc_matrix(n, n);
		b11 = alloc_matrix(n, n); b12 = alloc_matrix(n, n);
		b21 = alloc_matrix(n, n); b22 = alloc_matrix(n, n);
		c11 = alloc_matrix(n, n); c12 = alloc_matrix(n, n);
		c21 = alloc_matrix(n, n); c22 = alloc_matrix(n, n);
		m1 = alloc_matrix(n, n); m2 = alloc_matrix(n, n);
		m3 = alloc_matrix(n, n); m4 = alloc_matrix(n, n);
		m5 = alloc_matrix(n, n); m6 = alloc_matrix(n, n);
		m7 = alloc_matrix(n, n);
		aresult = alloc_matrix(n, n);
		bresult = alloc_matrix(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				a11[i][j] = a[i][j];
				a12[i][j] = a[i][j + n];
				a21[i][j] = a[i + n][j];
				a22[i][j] = a[i + n][j + n];
				b11[i][j] = b[i][j];
				b12[i][j] = b[i][j + n];
				b21[i][j] = b[i + n][j];
				b22[i][j] = b[i + n][j + n];
			}
		}
		plus(a11, a22, aresult, n);
		plus(b11, b22, bresult, n);
		strassen_multiply(aresult, bresult, m1, n);
		plus(a21, a22, aresult, n);
		strassen_multiply(aresult, b11, m2, n);
		minus(b12, b22, bresult, n);
		strassen_multiply(a11, bresult, m3, n);
		minus(b21, b11, bresult, n);
		strassen_multiply(a22, bresult, m4, n);
		plus(a11, a12, aresult, n);
		strassen_multiply(aresult, b22, m5, n);
		minus(a21, a11, aresult, n);
		plus(b11, b12, bresult, n);
		strassen_multiply(aresult, bresult, m6, n);
		minus(a12, a22, aresult, n);
		plus(b21, b22, bresult, n);
		strassen_multiply(aresult, bresult, m7, n);
		plus(m3, m5, c12, n);
		plus(m2, m4, c21, n);
		plus(m1, m4, aresult, n);
		plus(aresult, m7, bresult, n);
		minus(bresult, m5, c11, n);
		plus(m1, m3, aresult, n);
		plus(aresult, m6, bresult, n);
		minus(bresult, m2, c22, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				c[i][j] = c11[i][j];
				c[i][j + n] = c12[i][j];
				c[i + n][j] = c21[i][j];
				c[i + n][j + n] = c22[i][j];
			}
		}
		free_m(n, a11); free_m(n, a12);
		free_m(n, a21); free_m(n, a22);
		free_m(n, b11); free_m(n, b12);
		free_m(n, b21); free_m(n, b22);
		free_m(n, c11); free_m(n, c12);
		free_m(n, c21); free_m(n, c22);
		free_m(n, aresult); free_m(n, bresult);
		free_m(n, m1); free_m(n, m2);
		free_m(n, m3); free_m(n, m4);
		free_m(n, m5); free_m(n, m6);
		free_m(n, m7);
	}
}

//introducerea rezultatul inmultirii Strassen in array
void strassen(int ****matrices, int ***sizematrix, int *length)
{
	int index1, index2;
	scanf("%d%d", &index1, &index2);
	if (index1 >= *length || index2 >= *length || index1 < 0 || index2 < 0) {
		printf("No matrix with the given index\n");
		return;
	}
	int ***copy = *matrices;
	int **size_copy = *sizematrix;
	int m = size_copy[index1][0];
	int **c = alloc_matrix(m, m);
	strassen_multiply(copy[index1], copy[index2], c, m);

	copy[*length] = c;
	size_copy[*length] = (int *)malloc(sizeof(int) * 2);
	size_copy[*length][0] = m;
	size_copy[*length][1] = m;
	(*length)++;
}

int main(void)
{
	char query;
	int *capacity;
	int *length;
	int **sizematrix;
	int ***matrices;
	capacity = malloc(sizeof(int));
	length = malloc(sizeof(int));
	matrices = (int ***)malloc(sizeof(int **));
	sizematrix = (int **)malloc(sizeof(int *));
	*capacity = 0;
	*length = 0;
	while (1) {
		scanf(" %c", &query);
		if (query == 'Q') {
		    for (int i = 0; i < *length; i++) {
				for (int j = 0; j < sizematrix[i][0]; j++)
					free(matrices[i][j]);
				free(matrices[i]);
				free(sizematrix[i]);
			}
			free(matrices);
			free(sizematrix);
			free(capacity);
			free(length);
			break;
		} else if (query == 'L') {
			int extend = needs_extend(capacity, length);
			if (extend == 1) {
				matrices = extend_matrices(matrices, capacity);
				sizematrix = extend_size_matrix(sizematrix, capacity);
			}
			read_new_matrix(matrices, sizematrix, length);
		} else if (query == 'P') {
			print_matrix(matrices, sizematrix, length);
		} else if (query == 'C') {
			resize(&matrices, &sizematrix, length);
		} else if (query == 'D') {
			dimension(sizematrix, length);
		} else if (query == 'M') {
			int extend1 = needs_extend(capacity, length);
			if (extend1 == 1) {
				matrices = extend_matrices(matrices, capacity);
				sizematrix = extend_size_matrix(sizematrix, capacity);
			} mat_multiplication(&matrices, &sizematrix, length);
		} else if (query == 'T') {
			transpose(&matrices, &sizematrix, length);
		} else if (query == 'R') {
			logarithmic_pow(&matrices, &sizematrix, length);
		} else if (query == 'O') {
			matrix_sort(&matrices, &sizematrix, length);
		} else if (query == 'F') {
			free_matrix(&matrices, length, &sizematrix);
			int shrenk = needs_shrenk(capacity, length);
			if (shrenk == 1) {
				matrices = shrenk_matrix(matrices, capacity);
				sizematrix = shrenk_size_matrix(sizematrix, capacity);
			}
		} else if (query == 'S') {
			int extend2 = needs_extend(capacity, length);
			if (extend2 == 1) {
				matrices = extend_matrices(matrices, capacity);
				sizematrix = extend_size_matrix(sizematrix, capacity);
			} strassen(&matrices, &sizematrix, length);
		} else {
			printf("Unrecognized command\n");
		}
	}
	return 0;
}
