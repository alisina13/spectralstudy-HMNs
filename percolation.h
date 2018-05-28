#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "Random.h"

using namespace std;
/*
namespace Output {
    int realizations = 0;
    int normalization = 1;
    vector<int64_t> abscissa;
    vector<int64_t> data;
    vector<int64_t> csizes;
    void export_distribution(const char* filename)
    {
        size_t s = abscissa.size();
        ofstream f;
        f << scientific;
        f.open(filename);
        for(size_t i=0; i<s-1; ++i) {
            double interval = abscissa[i+1]-abscissa[i]; //implicit cast
            f << static_cast<double>(abscissa[i])/normalization << " ";
            f << static_cast<double>(data[i])/(interval*realizations) << "\n";
        }
        f.close();
    }
    void export_largest(const char* filename)
    {
        ofstream f;
        f << scientific;
        f.open(filename);
        for(auto& a : csizes){f<<static_cast<double>(a)/realizations<<"\n";}
        f.close();
    }
}
*/

template<typename T>
struct Node {
    vector<T> neighbors;
};



class Ensemble {
    Random ra;
    vector<Node<int>> nodelist;
    vector<Node<int>> newnl;
    vector< unsigned int > Giant_Component;
    vector<Node<int>> GC;
    vector<Node<int>> RelabledGC;
    int N;
    int steps;
    int fixed_links;
    public:

    inline double doub() {return ra.doub();}
    Ensemble(uint64_t seed, int N, int steps, int fl)
        : ra(Random(seed)), N(N), steps(steps), fixed_links(fl) {}
    void destroy_network()
    {
      nodelist.clear();
      //newnl.clear();
      //RelabledGC.clear();
    }
    void destroy_Percolated_GC_network()
    {
      newnl.clear();
      RelabledGC.clear();
    }
    void build_network_fixed()
    {
        /*Define network structure*/
        int n_of_blocks = (int)pow(2.0, (double)steps);
        int final_block_size = (int)(N / n_of_blocks);
        //	cout << n_of_blocks << " final blocks of size " << final_block_size << endl;
        /*Create nodes*/
        nodelist.reserve(N);					/*Allocate memory for array of node-objects*/
        for(int i=0; i<N; ++i) nodelist.emplace_back();
        /*Establish links*/
        int nlinks = 0;
        int corrections = 0;
        /*distance 0, inside blocks*/
        for(int block=0; block<n_of_blocks; block++) {
            int base = block*final_block_size;
            int added_links;
            added_links=0;
            for(int i=0; i<final_block_size; i++) {
                for(int j=i+1; j<final_block_size; j++) {
                    nodelist[i+base].neighbors.push_back((j+base));
                    nodelist[j+base].neighbors.push_back((i+base));
                    added_links++;
                }
            }
            nlinks+=added_links;
        }
        /*distances > 0, outside blocks*/
        int current_block_size = final_block_size;
        int current_n_of_blocks = n_of_blocks;
        for(int level=1; level<steps+1; level++) {
            for(int first=0; first<(int)(current_n_of_blocks/2); first++) {
                int established_links=0;
                while(established_links<fixed_links) {
                    int base1 = 2*first*current_block_size;
                    int base2 = base1+current_block_size;
                    double r1 = doub();
                    double r2 = doub();
                    int n1 = base1 + floor( r1*current_block_size );
                    int n2 = base2 + floor( r2*current_block_size );
                    bool already_linked=false;
                    for(int i=0; i<nodelist[n1].neighbors.size(); i++) {
                        if(nodelist[n1].neighbors[i] == n2)	{
                            already_linked=true;
                            break;
                        }
                    }
                    if(already_linked==false) {
                        nodelist[n1].neighbors.push_back(n2);
                        nodelist[n2].neighbors.push_back(n1);
                        established_links++;
                        nlinks++;
                    }
                }
            }
            current_block_size = (int)(current_block_size*2);
            current_n_of_blocks = (int)(current_n_of_blocks/2);
        }
        if(current_block_size!=N || current_n_of_blocks!=1) {
            cout << "Error building network" << endl;
            //exit(1);
        }
        //	cout << "Network generated: " << nlinks << " links"  << endl;
        //	cout << "corrections = " << corrections << endl;

    }

    vector<Node<int>> release_link_diluted_copy(double premove)
    {
        vector<int> fr; vector<int> to; vector<int> rm;
        for(int i=0; i<N; ++i) {
            for(auto& j : nodelist[i].neighbors) {
                if(i<j) {
                    fr.push_back(i); to.push_back(j); rm.push_back(0);
                }
            }
        }
        size_t listsize = rm.size();
        size_t nremove = listsize*(1.-premove);
        //cout << nremove << "!" << endl;
        if(nremove <= listsize) {
            size_t removed = 0;
            while(removed<nremove) {
                size_t target = floor(ra.doub()*listsize);
                if(!rm[target]) { //check that it is not marked for removal already
                    rm[target] = 1;
                    ++removed;
                }
            }
            std::cout<<removed<<" "<<"Number of links has been removed"<<std::endl;
        }
        //for(int i=0; i<listsize; ++i) {
        //	cout << fr[i] << " " << to[i] << " " << rm[i] << endl;
        //}
       newnl.reserve(N);
        for(int i=0; i<N; ++i) newnl.emplace_back();
        for(int i=0; i<listsize; ++i) {
            if(!rm[i]) {
                newnl[fr[i]].neighbors.push_back(to[i]);
                newnl[to[i]].neighbors.push_back(fr[i]);
            }
        }
        return newnl;
    }

    void Giant_Component_of_network()
    {
	     std::cout<<"start colouring different regions in network"<<std::endl;
	     unsigned int start_point = 0;
	     int maximum_colours = 0;
	     unsigned int number_of_colours(N);// The size of the percolated network should be given
       //vector<unsigned int> Giant_Component;
	     vector<int>coloured_region;
       coloured_region.resize(number_of_colours);
	     for(unsigned int n=0;n!=N;++n)
	        {
		          coloured_region[n] = -1;
	        }

		   for(unsigned int region_colour=0; region_colour!=number_of_colours; ++region_colour)
		      {
			         start_point = std::distance(coloured_region.begin(), std::find(coloured_region.begin(), coloured_region.end(), -1));
			         //nodelist[start_point].region=region_colour;
               if(start_point == N)
					          break;

			         vector<unsigned int> neighbours_list(1,start_point);
			         coloured_region[start_point] = region_colour;
			         while (!neighbours_list.empty())
			            {
				                vector<unsigned int> tmp_neighbours;
				                    for (auto x = neighbours_list.begin(); x != neighbours_list.end(); ++x)
				                        {
					                             for (auto y = newnl[*x].neighbors.begin(); y!= newnl[*x].neighbors.end(); ++y)
					                                  {
						                                        if (coloured_region[*y] == -1)
						                                              {
							                                                     tmp_neighbours.push_back(*y);
							                                                     coloured_region[*y] = region_colour;
							                                                     //nodelist[*y].region=region_colour;
						                                               }
					                                  }
				                        }
				                neighbours_list = tmp_neighbours;
			           }
            //std::cout<<"number of components: "<<region_colour<<std::endl;
		    }

		  std::cout<<"different regions in network has been coloured"<<std::endl;
		  std::cout<<"start finding giant component"<<std::endl;
		  maximum_colours = *max_element(coloured_region.begin(), coloured_region.end());
      std::cout<<maximum_colours<<std::endl;
		  vector<unsigned int> cnt_colours(maximum_colours+1,0);
		  for(unsigned int i = 0;i!=N;++i )
			{
				//std::cerr << "\t" <<  coloured_region[i] << "\n";
				cnt_colours.at(coloured_region[i])+=1;
			}
			int colour_of_Giant_Component = std::distance(cnt_colours.begin(),max_element(cnt_colours.begin(), cnt_colours.end()));
      //Giant_component.reserve(N);
		  for(unsigned int i = 0;i!=N;++i )
			{
				if(coloured_region[i]==colour_of_Giant_Component )
					{
					  Giant_Component.emplace_back(i);
					}
			}

		 std::cout<<"giant_component found"<<std::endl;
         std::cout<<"size of Gian Component is: "<<Giant_Component.size()<<std::endl;

     }

    void release_GC_of_network()
    {

        GC.reserve(Giant_Component.size());
        for(int i=0; i!=Giant_Component.size(); ++i)
        {
            GC.push_back(newnl[Giant_Component[i]]);
            GC[i].neighbors = newnl[Giant_Component[i]].neighbors;
        }
    }


    void relable_GC_of_network()
    {
      RelabledGC.reserve(Giant_Component.size());
      for(int i=0; i!=Giant_Component.size(); ++i)
      {RelabledGC.emplace_back();}

      for(int i=0; i!=Giant_Component.size(); ++i)
      {
          int degree=newnl[Giant_Component[i]].neighbors.size();
          for(int j=0; j!=degree; ++j)
          {
            RelabledGC[i].neighbors.push_back(std::distance( Giant_Component.begin() ,std::find(Giant_Component.begin(), Giant_Component.end(), newnl[Giant_Component[i]].neighbors[j]) ) ) ;
          }
      }
    }

    void degree_distribution() {

      vector<int> deg_dist;

      for (unsigned int i = 0; i != Giant_Component.size(); ++i) {
          deg_dist.push_back(RelabledGC[i].neighbors.size());
      }
      auto it = max_element(std::begin(deg_dist), std::end(deg_dist));
      std::cout<<"largest degree of GC of percolated network   " << *it <<std::endl;
    }




    vector<Node<int> >  & return_network_ref(){
      return nodelist;
    }

    vector<Node<int> >  & return_Percolated_network_ref(){
      return newnl;
    }

    vector<Node<int> >  & return_GC_of_Percolated_network_ref()
    {
      return GC;
    }

    vector<Node<int> >  & return_RelabledGC_of_Percolated_network_ref()
    {
      return RelabledGC;
    }


    void output_edgelist( vector<Node<int>> & n )
    {
        for(int i=0; i<N; ++i) {
            for(auto& j : n[i].neighbors) {
                if(i<j) cout << i << " " << j << endl;
            }
        }
    }

    void output_edgelist_GC( vector< Node<int> > & n )
    {
      std::ofstream file("GC.txt");
      for(int i=0; i<Giant_Component.size(); ++i) {
          //for(auto& j : n[i].neighbors) {
              /*if(i<j)*/ file << Giant_Component[i] <</* " " << j <<*/ endl;
        //  }
      }
    }


    void output_edgelist_RelabledGC( vector< Node<int> > & n )
    {
      for(int i=0; i<Giant_Component.size(); ++i) {
          for(auto& j : n[i].neighbors) {
              if(i<j) cout << i << " " << j << endl;
          }
      }
    }



};
