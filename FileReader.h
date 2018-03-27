#ifndef FILEREADER_HPP
#define FILEREADER_HPP

#include "Types.hpp"
#include <string>




//*******************************************************************************************************************
/*! Class for reading configuration from a file
*
* Configuration File Syntax:
*   - everything between a '#' character and the beginning of the next line is ignored
*   - lines can be empty
*   - line contain parameter key followed by white-spaces followed by their value
*
*  All possible keys (and the datatype of the value) have to be registered first:
*   - for example usage have a look at the FileReaderTest
*
*
*
*  This Skeleton is just a suggestion. If you are familiar with template programming
*  a possible alternative would be a version that does not need registering of all keys
*  and has a getter function like the following:
*      template<typename ExpectedType>
*      ExpectedType getValue( const std::string & key );
*/
//*******************************************************************************************************************
class FileReader
{
public:

	//register a new parameter with name key and initial int value
	void registerIntParameter( const std::string & key, int init = 0 );

	//register a new parameter with name key and initial double value
	void registerRealParameter( const std::string & key, real init = 0 );

	//register a new parameter with name key and initial string value
	void registerStringParameter( const std::string & key, const std::string & init = "" );

	//set a value for the key string with value in
	void setParameter( const std::string & key, const std::string & in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, real in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, int in );

	// get the int value of key 
	inline int getIntParameter( const std::string & key ) const;

	// get the double value of key 
	inline real getRealParameter( const std::string & key ) const;

	// get the string value of key 
	inline std::string getStringParameter( const std::string & key ) const;

	//try to read all registered parameters from file name
	bool readFile( const std::string & name );

	//print out all parameters to std:out
       	void printParameters() const;
	
	//checking wether the value is found
	//bool checkStrings(const std::string & key);
	//bool checkInt(const std::string & key);
	//bool checkReal(const std::string & key);
	
	
private:
	std::map<std::string,int>intmap;
	std::map<std::string,real>realmap;
	std::map<std::string,std::string> stringmap;

};




inline int FileReader::getIntParameter(const std::string &key) const
{
 //if( intmap.find(key)== intmap.end() )
   //   return "noKey" ; 
    //else
      return intmap.find(key)->second;
}

inline real FileReader::getRealParameter(const std::string &key) const
{
   //if( realmap.find(key)== realmap.end() )
     // return "noKey" ; 
    //else
      return (realmap.find(key))->second;
}

inline std::string FileReader::getStringParameter(const std::string &key) const
{
    if( stringmap.find(key)== stringmap.end() ) 
      return "noKey" ; 
    else
      return (stringmap.find(key))->second;
}



#endif //FILEREADER_HPP

