
#include "FileReader.h"
#include <map>
#include<fstream>
#include<iostream>
#include<sstream>
#include<algorithm>
using namespace std;


void FileReader::registerIntParameter(const std::string &key, int init)
{

     intmap.insert ( std::pair<const std::string,int>(key,init) );

}

void FileReader::registerRealParameter(const std::string &key, real init)
{

     realmap.insert ( std::pair<const std::string,real>(key,init) );

}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{

     stringmap.insert ( std::pair<const std::string,const std::string>(key,init) );

}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
    stringmap[key]=in;

}

void FileReader::setParameter(const std::string &key, real in)
{
    realmap[key]=in;
}

void FileReader::setParameter(const std::string &key, int in)
{
    intmap[key]=in;
}


bool FileReader::readFile(const std::string &name)
{
   ifstream infile(name.c_str());
   if (!infile) {
       std::cout<<"opening file failed"<<endl;
       exit(0);
       }
       string line;
       string key;
       int intvalue;
       real realvalue;
       string stringvalue;

       while(getline(infile,line)){
		// finding # and deleting strings after that
		char symbol('#');
		size_t found=line.find(symbol);
		if (found!=std::string::npos)
		line.erase(line.begin() + found,line.end());


	     // skipping whitespace
	        line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));


		stringstream	ss(line);
		ss>>key;

	       	for(sii it=intmap.begin(); it!=intmap.end(); ++it)
			 if (key.compare((*it).first)==0){
				                 	  ss>>intvalue;
						 	  intmap[key]=intvalue;
							 }
		for(sri it=realmap.begin(); it!=realmap.end(); ++it)
               		 if (key.compare((*it).first)==0){
		  				           ss>>realvalue;
					                   realmap[key]=realvalue;
						  	  }
		for(ssi it=stringmap.begin(); it!=stringmap.end(); ++it)
                	 if (key.compare((*it).first)==0){
		   					     ss>>stringvalue;
							     stringmap[key]=stringvalue;
							     }


                                   }

		   return true;

}



void FileReader::printParameters() const
{
 for(sii it=intmap.begin(); it!=intmap.end(); ++it)
 cout<<(*it).first<<" "<<(*it).second<<endl;

 for(sri it=realmap.begin(); it!=realmap.end(); ++it)
 cout<<(*it).first<<" "<<(*it).second<<endl;

 for(ssi it=stringmap.begin(); it!=stringmap.end(); ++it)
 cout<<(*it).first<<" "<<(*it).second<<endl;

}
