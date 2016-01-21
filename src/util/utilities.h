#ifndef _UTILITIES_H_
#define _UTILITIES_H_
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>
#define DOM_P -1
#define DOM_Q 1
#define DOM_INCOMPARABLE 0
using namespace std;

float** AllocateDoubleArray(const unsigned row, const unsigned col) {
  float** matrix = new float*[row];
  for (unsigned i = 0; i < row; i++)
    matrix[i] = new float[col];

  return matrix;
}

void FreeDoubleArray(const unsigned row, float** matrix) {
  for (unsigned i = 0; i < row; i++)
    delete[] matrix[i];
  delete[] matrix;
}

void PrintSkyline(const vector<int> &sky) {
//  printf(" ids:");
  for (uint32_t i = 0; i < sky.size(); ++i) {
    printf( " %d", sky[i] );
  }
  printf( "\n" );
}

bool CompareTwoLists(vector<int>& list1, vector<int>& list2, bool print_missing) {
  bool flag = true;
  if ( list1.size() == list2.size() ) {
    sort( list1.begin(), list1.end() );
    sort( list2.begin(), list2.end() );
    for (unsigned i = 0; i < list1.size(); i++) {
      if ( list1[i] != list2[i] ) {
        flag = false;
        break;
      }
    }
  } else {
    flag = false;
  }

  if (!flag && print_missing) {
    sort( list1.begin(), list1.end() );
    sort( list2.begin(), list2.end() );
    printf( "list1 missing:" );
    for (uint32_t i = 0; i < list2.size(); ++i) {
      if ( !std::binary_search( list1.begin(), list1.end(), list2[i] ) ) {
        printf( " %d", list2[i] );
      }
    }
    printf( "\n" );
    printf( "list2 missing:" );
    for (uint32_t i = 0; i < list1.size(); ++i) {
      if ( !std::binary_search( list2.begin(), list2.end(), list1[i] ) ) {
        printf( " %d", list1[i] );
      }
    }
    printf( "\n" );
  }

  return flag;
}

vector<string> &my_split(const string &s, char delim, vector<string> &elems) {
  stringstream ss( s );
  string item;
  while ( std::getline( ss, item, delim ) ) {
    elems.push_back( item );
  }
  return elems;
}

vector<string> my_split(const string &s, char delim) {
  vector<string> elems;
  my_split( s, delim, elems );
  return elems;
}

//reads a number of type T from a string.
template<class T> bool from_string(T& t, const std::string& s,
    std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss( s );
  return !(iss >> f >> t).fail();
}

//used to split each line of the input into strings and convert them to floats
vector<float> split(string ins, bool has_line_numbers) {
  stringstream ss( ins );
  string s;
  float f;
  vector<float> returnvals;
  bool ignore_first = false;
  if ( has_line_numbers ) {
    ignore_first = true;
  }
  while ( getline( ss, s, ',' ) ) {
    if ( ignore_first ) {
      ignore_first = false;
    } else {
      if ( from_string<float>( f, std::string( s ), std::dec ) ) {
        returnvals.push_back( f );
      } else {
        cout << "PARSE FAILED" << endl;
      }
    }
  }
  return returnvals;
}

vector<float> split_int(string ins, bool line_numbers) {
  stringstream ss( ins );
  string s;
  vector<float> returnvals;
  bool ignore_first = false;
  if ( line_numbers ) {
    ignore_first = true;
  }
  while ( getline( ss, s, ',' ) ) {
    if ( ignore_first ) {
      ignore_first = false;
    } else {

      int numb;
      istringstream( s ) >> numb;
      returnvals.push_back( (float) numb );

    }
  }
  return returnvals;
}

vector<vector<float> > split_data(vector<vector<float> > input, int d) {
  vector<vector<float> > ret;
  for (uint32_t i = 0; i < input.size(); i++) {
    vector<float> tempVec;
    for (int j = 0; j < d; j++) {
      tempVec.push_back( input.at( i ).at( j ) );
    }
    ret.push_back( tempVec );
  }
  return ret;
}

//reads the file at filename and optionally removes line numbers and normalizes
//the input to a range of [0,1]
vector<vector<float> > read_data(const char *filename, bool has_line_numbers,
    bool normalize) {
  ifstream file( filename );
  if ( !file.good() ) {
    printf( "Can't find '%s' file\n", filename );
    exit(EXIT_FAILURE);
  }
  string line;
  vector<vector<float> > lines;
  vector<float> tempVec;
  getline( file, line );
  tempVec = split( line, has_line_numbers );
  //how many entries do we have?
  int size = tempVec.size();
  float max[size];
  float min[size];
  for (unsigned int i = 0; i < tempVec.size(); i++) {
    max[i] = tempVec.at( i );
    min[i] = tempVec.at( i );
  }
  lines.push_back( tempVec );
  while ( getline( file, line ) ) {
    tempVec = split( line, has_line_numbers );
    for (unsigned int i = 0; i < tempVec.size(); i++) {
      if ( tempVec.at( i ) > max[i] )
        max[i] = tempVec.at( i );
      if ( tempVec.at( i ) < min[i] )
        min[i] = tempVec.at( i );
    }
    lines.push_back( tempVec );
  }
  if ( normalize ) {
    for (unsigned int i = 0; i < lines.size(); i++) {
      tempVec = lines.at( i );
      for (unsigned int j = 0; j < tempVec.size(); j++) {
        tempVec.at( j ) = (tempVec.at( j ) - min[j]) / (max[j] - min[j]);
      }
      lines.at( i ) = tempVec;
      //printf("%f",tempVec.at(0));
    }
  }
  return lines;
}

bool point_equal(vector<float> p_stop, vector<float> next) {
  for (uint32_t i = 0; i < p_stop.size(); i++) {
    if ( p_stop.at( i ) != next.at( i ) ) {
      return false;
    }
  }
  return true;
}

//used to determine the dominance in the sequential BNL algorithm
int dominates(vector<float> p, vector<float> q) {

  // Initialize dominator to incomparable
  // Note: This will be the result for two points that are equal
  // in each dimension
  int dominator = DOM_INCOMPARABLE;
  int dimensionality = p.size();
  // Go through dimensions to determine domination
  for (int i = 0; i < dimensionality; i++) {

    if ( p.at( i ) < q.at( i ) ) {
      if ( dominator == DOM_Q ) {
        // If q was set as dominating, but p is better
        // in a given dimension, then the points are incomparable.
        dominator = DOM_INCOMPARABLE;
        break;
      } else {
        // Otherwise just set P as dominating.
        dominator = DOM_P;
      }
    } else if ( p.at( i ) > q.at( i ) ) {
      // Same as above, just reversed
      if ( dominator == DOM_P ) {
        dominator = DOM_INCOMPARABLE;
        break;
      } else {
        dominator = DOM_Q;
      }
    }
  }

  return dominator;
}

//converts a multi vector to a single vector by putting the floats from
//the multi vector into one big vector
vector<float> to_single_vector(vector<vector<float> > dataset) {
  vector<float> tempVec;
  vector<float> result;
  int k = dataset.front().size();
  for (unsigned int i = 0; i < dataset.size(); i++) {
    tempVec = dataset.at( i );
    for (int j = 0; j < k; j++) {
      result.push_back( tempVec.at( j ) );
    }
  }
  return result;
}

//converts a multi vector to a matrix of floats
void redistribute_data(vector<vector<float> > datasetv, float** dataset) {
  for (unsigned int i = 0; i < datasetv.size(); i++) {
    float* x = dataset[i];
    vector<float> next = datasetv.at( i );
    for (unsigned int j = 0; j < datasetv.front().size(); j++) {
      x[j] = next[j];
    }
  }
}

#endif /* _UTILITIES_H_ */
