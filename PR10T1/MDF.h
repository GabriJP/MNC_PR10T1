/*


*/

#define CONTOUR_TYPE1 1
#define CONTOUR_TYPE2 2
#define CONTOUR_TYPE3 3
#define CONTOUR_TYPE4 4


class MDF{

public:
	MDF();
	~MDF();
	int solve(	const int N, const double h, 
				double *y, const double *params, 
				const int type, const double ca, const double cb);

private:
	void ensamble_kernel(const int N, const double h, const double *params, 
					double *LHS, double *RHS);

	// Contorno ya e yb
	void ensamble_contour1(const int N, const double h, double *LHS, double *RHS, double ya, double yb);

	// Contorno ya e y'b
	void ensamble_contour2(const int N, const double h, double *LHS, double *RHS, double ya, double y1b);

	// Contorno y'a e yb
	void ensamble_contour3(const int N, const double h, double *LHS, double *RHS, double y1a, double yb);

	// Contorno y'a e y'b
	void ensamble_contour4(const int N, const double h, double *LHS, double *RHS, double y1a, double y1b);

};