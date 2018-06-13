#ifdef __cplusplus
extern "C" void dgetrf_(
	const int &m,			// (input)
	const int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	int *ipiv,			// ipiv[min(m,n)] (output)
	int &info			// (output)
	);
#else /* ! __cplusplus */
void dgetrf_(
	const int *m,			/* (input) */
	const int *n,			/* (input) */
	double *a,			/* a[n][lda] (input/output) */
	const int *lda,			/* (input) */
	int *ipiv,			/* ipiv[min(m,n)] (output) */
	int *info			/* (output) */
	);
#endif /* ! __cplusplus */

