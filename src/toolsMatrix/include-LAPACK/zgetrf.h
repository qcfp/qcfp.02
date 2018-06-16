#ifdef __cplusplus
extern "C" void zgetrf_(
	const int &m,			// (input)
	const int &n,			// (input)
	complex<double> *a,		// a[n][lda] (input/output)
	const int &lda,			// (input)
	int *ipiv,			// ipiv[min(m,n)] (output)
	int &info			// (output)
	);
#else /* ! __cplusplus */
void zgetrf_(
	const int *m,			/* (input) */
	const int *n,			/* (input) */
	complex_double *a,		/* a[n][lda] (input/output) */
	const int *lda,			/* (input) */
	int *ipiv,			/* ipiv[min(m,n)] (output) */
	int *info			/* (output) */
	);
#endif /* ! __cplusplus */

