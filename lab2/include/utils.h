#include <fstream>
#include <vector>
#include <iostream>

#ifndef LAB2_MATRIX_H
#define LAB2_MATRIX_H

#endif //LAB2_MATRIX_H

using namespace std;

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

void generateRandomMatrix(const string& fname, int n){
    ofstream file(fname);
    double val;

    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            val = 1.0*rand() / RAND_MAX;
            file << val << " ";
        }
        file << "\n";
    }

    file.close();
}

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