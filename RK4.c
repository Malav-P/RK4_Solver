#include<stdio.h>
#include<stdlib.h>
#include <math.h>


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


/* iterate is a function to iterate the state vector, y_n+1 = y_n + (h/6)(a + 2*b + 2*c +d), in other words the RK order 4 algorithm

---INPUTS -------------------------------------------------------------------------------------

t : a double representing the time variable in the sytem
h : a double representing the time step between iterations
x : a pointer variable pointing to the first element of the state vector containing elements of type double
size: an int variable representing the number of elements in the state vector
f : a pointer variable to the function of the RHS of the equation xdot = func(t, x)

---OUTPUTS-------------------------------------------------------------------------------------

x_iter : a pointer variable to the array that describes the iterated state vector x_n+1 = x_n + (h/6)(a + 2*b + 2*c +d)

---BEGIN FUNCTION DECLARATION AND BODY----------------------------------------------------------------------------*/

double* iterate(double t, double h, double* x, int size, double* (*f)(double, double*, int)){
	int length = size;

	double* a = (*f)(t, x, size);
    
	double temp[length];
	double temp1[length];
	double temp2[length];
	int i;

	for (i = 0; i<length; ++i){
		temp[i] = *(x+i) + (h/2)*(*(a+i));
	}
	double* b = (*f)(t + h/2, temp, size);
	
	
	for (i = 0; i<length; ++i){
		temp1[i] = *(x+i) + (h/2)*(*(b+i));
	}
	double* c = (*f)(t + h/2, temp1, size);

	
	for (i = 0; i<length; ++i){
		temp2[i] = *(x+i) + (h)*(*(c+i));
	}
	double* d = (*f)(t + h, temp2, size);

    double* x_iter;
	x_iter = (double*)calloc(size, sizeof(double));

	for (i = 0; i<length; ++i){
		x_iter[i] = *(x+i) + (h/6)*(*(a+i)+2**(b+i)+2**(c+i)+*(d+i));
	}
    free(a);
    free(b);
    free(c);
    free(d);
	return x_iter;


}

/*---END FUNCTION DECLARATION AND BODY----------------------------------------------------------------------------*/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/* solve_ivp is a function to solve an IVP using RK4 method

---INPUTS -------------------------------------------------------------------------------------

x_0 : an array of doubles representing the initial state vector
size: an int representing the number of elements in the state vector
tspan: an array of two elements of type double that represents the time interval of integration
h : a double representing the time step between iterations
f : a pointer variable to the function of the RHS of the equation xdot = func(t, x)

---OUTPUTS-------------------------------------------------------------------------------------

answer : a pointer to an array of pointers that describe the matrix of solutions. Each row of answer is the state vector at time t

---BEGIN FUNCTION DECLARATION AND BODY-----------------------------------------------------------------------------*/

double** solve_ivp(double x_0[], int size, double t_span[], double h, double* (*f)(double, double*, int)){
    int no_of_points = ceil((t_span[1]-t_span[0])/h);
    
    double** answer;
    answer = (double**)calloc(no_of_points, sizeof(double*));
    int i;
    for(i=0;i<no_of_points;++i){
        answer[i] = (double*)calloc(size, sizeof(double));
    }

	int j;
	int k;
    
    for (j=0;j<size;j++){
        answer[0][j] = x_0[j];
    }

    double temp_var = t_span[0];
	for ( k = 1; k < no_of_points; ++k ){
		double* temp = iterate(temp_var, h, answer[k-1], size, &(*f));
		printf("state vector at t = %f: {", temp_var + h);
		for (j = 0; j<(size-1); ++j){
		answer[k][j] = *(temp + j);
		printf("%f, ", answer[k][j]);
		}
		answer[k][size-1] = *(temp + (size-1));
		printf("%f}\n", answer[k][size-1]);
		temp_var += h;
		free(temp);
	}
	return answer;
}

/**************************************************************************************************************************************/

/* dxdt is a function to return the right-hand-side of the equation xdot = f(t, x)

---INPUTS -------------------------------------------------------------------------------------

t : a double representing the time variable in the sytem
x : a pointer variable pointing to the first element of the state vector containing elements of type double
size: an int variable representing the number of elements in the state vector

---OUTPUTS-------------------------------------------------------------------------------------

f : a pointer variable to the array that describes the RHS of xdot = f(t, x)

-----------------------------------------------------------------------------------------------*/

double* dxdt(double t, double* x, int size){               /* double* states that dxdt returns a pointer that points to a double type*/

	double* f;                                             /* f is a pointer to &f[0], the first element in the array of doubles.*/
	f = (double*)calloc(size, sizeof(double));	           /* allocating space in memory for the array. f points to the first cell in the array*/
	f[0] = *(x+1);
	f[1] = -2*0.1*10**(x+1) - 100**(x);
	f[2] = 3*t*t;
	f[3] = 2*t;

	return f;                                              /* return the pointer that points to the first element in the array of doubles*/
}

/*****************************************************************************************************************************************/


/*---END FUNCTION DECLARATION AND BODY-----------------------------------------------------------------------------*/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


int main(){
    double y_0[] = {0, 1, 0, 0};                                /* the state vector of initial conditions*/
    double t_s[2] = {0, 2};                                     /* time span of integration*/
    double h = 0.0001;                                          /* time step of integrator*/


    int siz = sizeof(y_0)/sizeof(double);                  /* length of state vector*/
    double** ans = solve_ivp(y_0, siz, t_s, h, &dxdt);

    printf("program complete, press enter to exit\n");
    getchar();

	int no_of_points = ceil((t_s[1]-t_s[0])/h);
    for(int i=0;i<no_of_points;++i){
        free(ans[i]);
    }
	free(ans);
    
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/