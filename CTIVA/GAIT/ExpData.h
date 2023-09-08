#ifndef EXPDATA_H_
#define EXPDATA_H_

#include <vector>
#include <string>

using namespace std;


typedef struct _expdata_t {
	string name;
	vector<double> exp;
} expdata_t;

class ExpData
{
public:
	ExpData();
	~ExpData();


	vector<expdata_t>* getExpData(){return &mEXP;}
	bool getExpDataFromFile(char *filename);


private:
	int mSize;
	vector<expdata_t> mEXP;

};





#endif /* EXPDATA_H_ */
