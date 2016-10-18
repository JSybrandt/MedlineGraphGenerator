#include"Dict.h"

using namespace std;
Dict::Dict(string file) {
	ifstream dataFile(file);
	string line;
	bool firstLine = true;
	while (getline(dataFile,line)) {
		if (firstLine) { firstLine = false; continue;}
		if (line.length() > 0) {
			stringstream s;
			s << line;
			string word;
			s >> word;
			vector<float> vectorData;
			float d;
			while (s >> d) {
				vectorData.push_back(d);
			}
			data[word] = Vec(vectorData);
		}
	}
	dataFile.close();
}

const Vec Dict::getVec(string word) {
	if (data.find(word) == data.end())
		return Vec();
	else
		return data[word];
}

const bool Dict::contains(string word){
    return data.find(word) != data.end();
}