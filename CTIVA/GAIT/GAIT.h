#ifndef GAIT_H_
#define GAIT_H_

#include <vector>
#include <string>
#include "CenData.h"
#include "Region.h"
#include "ExpData.h"
#include "Region.h"
#include "OPT.h"
#include "PHI.h"

using namespace std;

typedef struct _tgen_t {
	vector<double>T1Gen;
	vector<double>T2Gen;
} tgen_t;


class GAIT
{
public:
	GAIT();
	~GAIT();

	void Tgen(int NDIV,int N_sample,CenData* Data,Density* Den, bool Bound);
	void gait(ExpData* ExpData);
	vector<tgen_t>* getTgen(){return &mTGen;}

private:
	int N_tgen;
	vector<tgen_t> mTGen;

};



#endif /* GAIT_H_ */
