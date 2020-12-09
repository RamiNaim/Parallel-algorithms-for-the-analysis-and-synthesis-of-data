#include <fstream>
#include <vector>
#include <iostream>

#ifndef LAB2_MATRIX_H
#define LAB2_MATRIX_H

#endif //LAB2_MATRIX_H

using namespace std;

void readMatrix(const string& fname, vector<vector<double>>* buffer){
    ifstream file(fname, ios::in);

    string delimiter = " ";
    string literal;
    vector<double> row;

    string line;
    size_t pos;
    if (file.is_open()){
        while (getline(file, line)){
            while ((pos = line.find(delimiter)) != string::npos) {
                literal = line.substr(0, pos);
                row.push_back(  stod(literal) );
                line.erase(0, pos + delimiter.length());
            }
            buffer -> push_back(row);
            row.clear();
        }
    }

    file.close();

}