#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <time.h>
#include <iomanip>
#include "mkl.h"
#include "Random.h"
#include "percolation.h"
#include "NetworkCSR.h"
#include "FileReader.h"


using namespace std;


template<typename T>
void get_from_command_line(int argc, char **argv, string name, T &a) {
    string value;
    for(int i=1; i<argc; ++i) {
        string argument(argv[i]);
        size_t found = argument.find(name);
        if(found!=string::npos) {
            size_t eqsgn = argument.find("=");
            value = string(argument,eqsgn+1,argument.size());
            break;
        }
    }
    stringstream vss(value);
    vss >> a;
}

int main(int argc, char **argv)
{
    FileReader   conf;
    /* Declaration of FEAST variables */
    char         UPLO = 'F'; /* Type of matrix: (F means full matrix, L/U - lower/upper triangular part of matrix) */
    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

    /* Declaration of local variables */
    conf.registerRealParameter("prob");
    conf.registerIntParameter("seed");
    conf.registerIntParameter("N");
    conf.registerIntParameter("realizations");
    conf.registerIntParameter("steps");
    conf.registerIntParameter("fixed_links");
    conf.registerRealParameter("Emin");
    conf.registerRealParameter("Emax");
    conf.registerIntParameter("M0");
    conf.registerIntParameter("M");
    conf.registerIntParameter("loop");
    conf.registerIntParameter("info");
    conf.registerRealParameter("epsout");

   string inputfilepath = "FileReaderTestInput.txt";
   std::cout<<argv[0]<<argv[1]<<std::endl;
    // parseparameters
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: you can give an input file as parameter " << std::endl;

        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */


    }
    else {
    inputfilepath = argv[1];
  }



    //reading file
    conf.readFile( inputfilepath);
    conf.printParameters();

    double prob = conf.getRealParameter("prob");
    int seed = conf.getIntParameter("seed");
    int N = conf.getIntParameter("N");
    int steps = conf.getIntParameter("steps");
    int fixed_links = conf.getIntParameter("fixed_links");
    int realizations = conf.getIntParameter("realizations");
    double Emin = conf.getRealParameter("Emin"); /* Lower/upper bound of search interval [Emin,Emax] */
    double Emax = conf.getRealParameter("Emax");
    MKL_INT loop = conf.getIntParameter("loop");    /* Number of refinement loop */
    MKL_INT M0 = conf.getIntParameter("M0");   /* Initial guess for subspace dimension to be used */
    MKL_INT M = conf.getIntParameter("M");    /* Total number of eigenvalues found in the interval */
    double epsout = conf.getRealParameter("epsout");  /* Relative error on the trace */
    MKL_INT info = conf.getIntParameter("info");          /* Errors */


    Ensemble one(seed,N,steps,fixed_links); //(seed,N,steps,alpha)
    for(int i=0; i<realizations; ++i) {
        if(!(i%100)) cout << i << endl;
        cout << "Build AdjList" <<" with"<< " "<<N<<" "<<"nodes"<<" "<<"and"<<" "<<steps<<" "<<"steps"<< endl;
        one.build_network_fixed();
        //std::cout<<"output Network"<<std::endl;
        //one.output_edgelist(one.return_network_ref());

        one.release_link_diluted_copy(prob);           //percolation starts
        //std::cout<<"output Percolated Network"<<std::endl;
        //one.output_edgelist(one.return_Percolated_network_ref());
        one.destroy_network();
        one.Giant_Component_of_network();
        one.release_GC_of_network();
        std::cout<<"output GC of Percolated Network "<<std::endl;
        one.output_edgelist_GC(one.return_GC_of_Percolated_network_ref());
        one.relable_GC_of_network();
        one.degree_distribution();
        //std::cout<<"output Relabled GC of Percolated Network "<<std::endl;
        //one.output_edgelist_RelabledGC(one.return_RelabledGC_of_Percolated_network_ref());


        cout << "Build CSRMatrix" << endl;
        NetworkCSR<int,double> ncsr(one.return_RelabledGC_of_Percolated_network_ref(), 0);

        int newN= one.return_RelabledGC_of_Percolated_network_ref().size();
        double    *   E =  (double*)malloc(M*sizeof(double)) ;      /* Eigenvalues */
        double    *   X =  (double*)malloc(M*newN*sizeof(double));       /* Eigenvectors */
        double    *  res = (double*)malloc(M*sizeof(double)) ;      /* Residual */

        for (int i=0; i<M*one.return_RelabledGC_of_Percolated_network_ref().size(); i++)
        { X[i] = 0.;}

        for (int i=0; i<M; i++)
        {E[i] = 0.;}

        for (int i=0; i<M; i++)
        {res[i] = 0.;}

        //number of nonzeroes in upper triangular matrix, always including diagonal
        MKL_INT nnz = ncsr.get_n_of_links() + newN ;//The adjacency matrix is symmetric but the Normallized laplacian is not
        double  *a  = (double*)mkl_malloc(nnz*sizeof(double),64);
        MKL_INT *ja = (MKL_INT*)mkl_malloc(nnz*sizeof(MKL_INT),64);
        MKL_INT *ia = (MKL_INT*)mkl_malloc((newN+1)*sizeof(MKL_INT),64);

        ncsr.fill_csr3_matrix(ia,ja,a,nnz,one.return_RelabledGC_of_Percolated_network_ref());
      //  std::cout<<ia[0]<<" "<<ia[1]<<" "<<ia[2]<<" "<<ia[3]<<" "<<ia[4]<<" "<<ia[5]<<std::endl;
     //    std::cout<<ja[0]<<" "<<ja[1]<<" "<<ja[2]<<" "<<ja[3]<<" "<<ja[4]<<" "<<ja[5]<<" "<<ja[6]<<" "<<ja[7]<<" "<<ja[8]<<std::endl;
      //  std::cout<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<" "<<a[4]<<" "<<a[5]<<" "<<a[6]<<" "<<a[7]<<" "<<a[8]<<std::endl;
        //now invoke eigensolver with arguments n,ia,ja,a etc.
        /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
        cout << "Solve Standard Eigenvalue Problem" << endl;
        feastinit(fpm);

        fpm[0] =  1; /* Extended Eigensolver routines print runtime status to the screen. */



        fpm[64+0] = 1; /* No solver default */
        fpm[64+1] = 0;         /* modified -- minimum degree algorithm */
        fpm[64+3] = 0;         /* No iterative-direct algorithm */
        fpm[64+4] = 0;         /* No user fill-in reducing permutation */
        fpm[64+7] = 2;         /* Max numbers of iterative refinement steps */
        fpm[64+9] = 13;        /* Perturb the pivot elements with 1E-13 */
        fpm[64+10] = 0;        /* Do no use nonsymmetric permutation and scaling MPS */
        fpm[64+12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
        fpm[64+13] = 0;        /* Output: Number of perturbed pivots */
        fpm[64+17] = 0;       /* Input/Output: no Number of nonzeros in the factor LU */
        fpm[64+18] = 0;       /* Input/Output: no Mflops for LU factorization */
        fpm[64+19] = 0;        /* Output: Numbers of CG Iterations */
        fpm[64+34] = 0;        /* PARDISO use Fortran-style indexing for ia and ja arrays */


        /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
        dfeast_scsrev (&UPLO, &newN, a, ia, ja, fpm, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);

        mkl_free(a);
        mkl_free(ja);
        mkl_free(ia);
        printf("FEAST OUTPUT INFO %d \n",info);
        if ( info != 0 )
        {
          printf("Routine dfeast_scsrev returns code of ERROR: %i", (int)info);
          return 1;
        }

        //writing eigenvalues for network
        double threshold[] = {1e-3,2e-3,3e-3,1e-2,2e-2,3e-2,1e-1,2e-1,3e-1};
        vector<double> thresholds (threshold, threshold + sizeof(threshold) / sizeof(double) );
        ofstream processed_data;
        vector<double> counter (thresholds.size(),0);

        processed_data.open("processed_data.dat");
        /*
        processed_data<<"eigenvalues"<<" "<<"threshold1"<<" "<<"counters1"<<" "<<"threshold2"<<" "<<"counters2"<<" "<<"threshold3"<<" "<<"counters3"<<" "
        <<"threshold4"<<" "<<"counters4"<<" "<<"threshold5"<<" "<<"counters5"<<" "<<"threshold6"<<" "<<"counters6"<<" "<<"threshold7"<<" "<<"counters7"
        <<" "<<"threshold8"<<" "<<"counters8"<<" "<<"threshold9"<<" "<<"counters9"<<endl;*/

          for (int i=0; i<M; i++)
          {
            processed_data <<setprecision(8)<< E[i]<<" " ;
            for(int t = 0 ; t!=thresholds.size(); ++t)
            {
              for (int j=0; j<newN; j++)
              {
                if( abs( X[(i*newN)+j] ) > thresholds[t] )
                counter[t]+=1;
              }
            processed_data << counter[t] <<" "<< thresholds[t]<<" " ;
            counter[t] = 0;
           }
           processed_data <<endl;

         }
        processed_data.close();

        one.destroy_Percolated_GC_network();


    }


}
