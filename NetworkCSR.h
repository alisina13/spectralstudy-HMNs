#ifndef __NetworkCSR_H__
#define __NetworkCSR_H__
#include <iostream>
#include <vector>
using namespace std;

template<typename NodeInd, typename NodeVal>
class NetworkCSR {
    NodeInd n_of_nodes;
    NodeInd n_of_links;

public:
    NetworkCSR( vector<Node<NodeInd> > & n,
               NodeVal default_val) {
        n_of_nodes = n.size();
        n_of_links = 0;
        for(auto& a : n) {
            n_of_links += a.neighbors.size();
        }
        //std::cout<<n_of_nodes<<"------------"<<n_of_links<<std::endl;
      }

  void fill_csr3_matrix(NodeInd* ia, NodeInd* ja, NodeVal * a, NodeInd nnz,
  					   vector<Node<NodeInd> >& nodelist) {
     	//IndexT is the type use for indices (integer)
     	//ValueT for values (floating point)
      NodeInd N = nodelist.size();
      //std::cout<<"---------------------"<<N<<std::endl;
      NodeInd rowpt = 0;
      NodeInd valpt = 0;
      ia[rowpt] = 1;
      for(NodeInd row_id=0; row_id<N; ++row_id) {
          vector<int> &neighbors = nodelist[row_id].neighbors;
          int degree = neighbors.size();
          vector<NodeInd> els;
          //always store diagonal elements, even if =0, MKL likes it this way
          els.push_back(row_id);
          for(auto& nn: neighbors) {
              if(nn>row_id) els.push_back(nn);
          }
          //sort column indices is ascending order, MKL likes it this way
          sort(els.begin(),els.end());
          for(auto e: els) {
              ja[valpt] = e+1;//changed
              //cout<<row_id<<" ______"<<e<<endl;
              //for Laplacian matrix:
              a[valpt] = (e!=row_id)? -1.0 : static_cast<NodeVal>(degree);
              //for adjacency matrix:
              //a[valpt] = (e!=row_id)? 1.0 : 0.0;
              ++valpt;
              //cout<<valpt<<endl;
          }
          ia[++rowpt] = valpt+1;//changed
      }

  }


    NodeInd const & get_n_of_links(){
      return n_of_links;
    }

    NodeInd const & get_n_of_nodes(){
      return n_of_nodes;
    }

};
#endif
