#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <algorithm>
#include <omp.h>


using namespace std;


typedef unsigned char uchar;
using vvchar = vector<vector<char> >;
using pint = pair<int,int>;

vector<ofstream*> of_vec;
int range_start = 1;
string chr_num;
int compress_flag = 0;
int no_gen_map = 1;
int ignore_singletons = 0;
int ac_for_list = 100;
int matrix_flag = 2; //0 means run original mode, 1 means run only matrix mode, 2 means both
int subsample_pairs = 10; //For k greater than this value, 5*k pairs will be subsampled and other matrix values replaced with -2



vector<string> pullNames(string);
uchar ** convertVecToData(vvchar &);
double getGenCol(string line, int col, double & phys);

bool nearZero(double a) {
	if(a <= 0.000001 && a >= -0.000001) { return 1; }
	return 0;
}

vector<string> split(string s, string delim) {
	vector<string> out;
	if(s.size() == 0) { return out; }
	int pos = 0;
	string token;
	while ((pos = s.find(delim)) != string::npos) {
		token = s.substr(0,pos);
		out.push_back(token);
		s.erase(0,pos+delim.length());
	}
	out.push_back(s);
	return out;
}

void addToVector(vector<string> & a, string & b) {
	for(int i = 0; i < a.size(); i++) {
		if(a[i].compare(b) == 0) { return; }
	}
	a.push_back(b);

}

string getCategory(string & s) {
	if(s.find("HIGH") != string::npos) {
		return "HIGH";
	} else if(s.find("MODERATE") != string::npos) {
		return "MODERATE";
	} else if(s.find("LOW") != string::npos) {
		return "LOW";
	}
	return "MODIFIER";
}

string typeLabel(vector<string> & info_s) {
	string out;
	vector<string> vout;
	string category;
	int ann_found = 0;
	for(int i = 0; i < info_s.size(); i++) {
		vector<string> s = split(info_s[i],"=");
		if(s[0].compare("ANN") == 0) {
			category = getCategory(s[1]);
			vector<string> anns = split(s[1],",");
			for(int j = 0; j < anns.size(); j++) {
				vector<string> al = split(anns[j],"|");
				addToVector(vout,al[1]);
			}
			ann_found = 1;
		}
	}
	if(ann_found == 0) { return ""; }
	out += (category + ",");
	for(int i = 0; i < vout.size(); i++) {
		out += vout[i];
		if(i != vout.size() - 1) { out += ","; }
	}
	return out;


}

bool isCpg(vector<string> & info_s) {
	for(int i = 0; i < info_s.size(); i++) {
		if(info_s[i].compare("CPG_TAG") == 0) { return 1; }
	}
	return 0;

}

bool isIndel(vector<string> & info_s) {
	for(int i = 0; i < info_s.size(); i++) {
		if(info_s[i].compare("INDEL") == 0) { return 1; }
	}
	return 0;
}





class maxhapdata {
	public:
	vector<int> cutoffs;
	int maxhap, remaining;
	maxhapdata() {

	};
	maxhapdata(int size, list<int> & idl) {
		vector<int> t(size,-2);
		cutoffs = t;
		maxhap = -1;
		remaining = idl.size();
		for (auto i = idl.begin(); i != idl.end(); i++) {
			cutoffs[*i] = -1;
		}
	}
	
};

class hapset {
	public:
	int hap_id; 
	int snp_id;
	list<int> hap_idx_list;
	maxhapdata data_left;
	maxhapdata data_right;
	hapset(int hap, int snp, int hap_count, list<int>  hil) {
		hap_id = hap;
		snp_id = snp;
		hap_idx_list = hil;
		maxhapdata t1(hap_count,hap_idx_list);
		maxhapdata t2(hap_count,hap_idx_list);
		data_left = t1;
		data_right = t2;
	}
	hapset(int hap, int snp, int hap_count, list<int>  hil, int extra) {
		hap_id = hap;
		snp_id = snp;
		hap_idx_list = hil;
		hap_idx_list.push_back(extra);
		hap_idx_list.sort();
		maxhapdata t1(hap_count,hap_idx_list);
		maxhapdata t2(hap_count,hap_idx_list);
		data_left = t1;
		data_right = t2;
	}
	void fillLeft(int mn) {
		for(auto i = hap_idx_list.begin(); i != hap_idx_list.end(); i++) {
			if(data_left.cutoffs[*i] == -1) {
				data_left.cutoffs[*i] = mn;
			}
		}
	}
	void fillRight(int mx) {
		for(auto i = hap_idx_list.begin(); i != hap_idx_list.end(); i++) {
			if(data_right.cutoffs[*i] == -1) {
				data_right.cutoffs[*i] = mx;
			}
		}
	}
	void print() {
		cout << "Hap id: " << hap_id << endl;
		cout << "Snp id: " << snp_id << endl;
		cout << "Max: " << data_left.maxhap << " " << data_right.maxhap << endl;
		list<int>::iterator itr = hap_idx_list.begin();
		for(int i = 0; i < 10; i++) {
			cout << *itr << " " << data_left.cutoffs[i] << " " << data_right.cutoffs[i] << endl;
			++itr;
		}
	}

};

class snp_data {
	public: 
	int allele_count;
	int position;
	int phased;
	double genetic_position;
	string annotation;
	bool rare = 0;
	list<int> hap_ids;
	snp_data(int ac, int pos, int phs, list<int> hi) {
		allele_count = ac;
		position = pos;
		phased = phs;
		if (ac <= ac_for_list) { 
			rare = 1; 
			hap_ids = hi;
		}
	}
	snp_data(int ac, int pos, int phs, string ann, list<int> hi) {
		allele_count = ac;
		position = pos;
		phased = phs;
		annotation = ann;
		if (ac <= ac_for_list) { 
			rare = 1; 
			hap_ids = hi;
		}
	}
	bool isSingleton() {
		if(allele_count == 1) {
			return 1;
		}
		return 0;
	}
	bool isValidForComp() {
		if(allele_count >= 2 || (ignore_singletons == 0 && allele_count == 1)) {
			return 1;
		}
		//Add check for CpG flag
		return 0;
	}
};

void addGeno(vector<char> & hap, char hh, int idx) {
	uchar h = (unsigned char)((hh == '0') ? 0 : 1);
	if(compress_flag == 0) {
		hap.push_back(h);
		return;
	}
	int res = idx%8;
	int offset = pow(2,7-res);
	uchar add = (unsigned char)(int(h)*offset);
	if(res == 0) {
		hap.push_back(add);
	} else {
		hap[hap.size()-1] += add;
	}
	return;
}

class vcf_data {
	public:
	uchar ** data;
	int hap_count;
	int snp_count;
	int transpose;
	int compress;
	int gendist_minmax[2];
	vector<snp_data> snps;
	vector<string> indiv_list;
	char getHap(int hap, int snp) {
		if(compress_flag == 0) {
			return char(data[hap][snp]);
		}
		int idx = snp/8;
		int reg = pow(2,7-(snp%8));
		int ti = int(data[hap][idx])%(2*reg);
		ti /= reg;
		uchar out = (unsigned char)ti;
		//cout << hap << " " << snp << " " << idx << " " << reg << " " << ti << " " << int(out) << endl;
		//out /= reg;
		return char(out);
	}
	void readVcf(string filename) {
		istream *inf;
		ifstream infile;
		if(filename.size() == 0) { inf = &cin; }
		else { 
			infile.open(filename.c_str());
			inf = &infile;
		}
		//infile.open(filename.c_str());
		string line;
		vvchar hap_hold;
		bool vector_set = 0;
		int snp_num = 0;
		while(getline(*inf,line)) {
			if(line[0] == '#') {
				if(line[1] != '#') {
					indiv_list = pullNames(line);
				}
				continue;
			}
			if(line.size() == 0) {
				break;
			}
			stringstream s(line);
			string junk, ref, alt, hap, info;
			int phys;
			list<int> hap_ids;
			s >> chr_num >> phys >> junk >> ref >> alt >> junk >> junk >> info;
			vector<string> info_s = split(info,";");
			if(ref.size() != 1 || isIndel(info_s) || isCpg(info_s)) {
				continue;
			}
			int ac = 0;
			int phased = -1;
			int ind_idx = 0;
			//int snp_num = 0;
			string annotation = "";
			if(vector_set == 0) { cout << "VECTOR_SET: " << info_s.size() << endl; }
			if (info_s.size() != 0) {
				annotation = typeLabel(info_s);
			}
			s >> junk;
			while(s >> hap) {
				if(vector_set == 0) {
					vector<char> t;
					hap_hold.push_back(t);
					vector<char> t2;
					hap_hold.push_back(t2);
				}
				if(phased == -1) { phased = (hap[1] == '|') ? 1 : 0; }
				if(hap[0] != '0') { 
					ac++;
					hap_ids.push_back(ind_idx);
				}
				if(hap[2] != '0') {
					ac++;
					hap_ids.push_back(ind_idx+1);
				}
				addGeno(hap_hold[ind_idx],hap[0],snp_num);
				//ind_idx++;
				addGeno(hap_hold[ind_idx+1],hap[2],snp_num);
				ind_idx+=2;
			}
			snp_data snp(ac,phys,phased,annotation,hap_ids);
			snps.push_back(snp);
			vector_set = 1;
			snp_num++;
		}
		snp_count = snps.size();
		hap_count = hap_hold.size();
		data = convertVecToData(hap_hold);
	}
	void print() {
		cout << "Hap count: " << hap_count << endl;
		cout << "Snp count: " << snp_count << endl;
		for(int i = 0; i < snp_count; i++) {
			cout << "Snp " << i << ", position " << snps[i].position << ", genetic " << snps[i].genetic_position << " " << gendist_minmax[0] << " " << gendist_minmax[1] << " " ;
			for(int j = 0; j < hap_count; j++) {
				//cout << int(data[j][i]);
				cout << int(getHap(j,i));
			}
			cout << endl;
		}
	}
	vector<int> getRareIdx(int mn, int mx, int pos) {
		vector<int> out;
		for(int i = 0; i < snps.size(); i++) {
			if(snps[i].allele_count >= mn && snps[i].allele_count <= mx && snps[i].position >= pos) {
				out.push_back(i);
				//cout << i << endl;
			}
		}
		return out;
	}
	int getSnpDist(int base, int end) {
		return snps[base].position - snps[end].position;
	}
	double getSnpGm(int base, int end) {
		return snps[base].genetic_position - snps[end].genetic_position;
	}
	void addGenData(string filename, int offset) {
		ifstream gm;
		gm.open(filename.c_str());
		string line;
		vector<double> g0,g1;
		double a,b;
		getline(gm,line);
		while(gm) {
			getline(gm,line);
			b = getGenCol(line,offset,a);
			if(gm.eof()) { break; }
			if(b <= -0.5) {
				cerr << "Error with genetic map\n";
				cerr << a << " " << b << " " << line << endl;
				exit(0);
			}
			//if(g1.size() == 0 || b != g1[g1.size()-1]) {
			if(!nearZero(b) && (g1.size() < 2 || b != g1[g1.size()-1])) {
				g0.push_back(a);
				g1.push_back(b);
				//cout << a << " " << b << endl;
			}
		}
		cerr << "Finding min\n";
		for(int ii = 1; ii < g0.size(); ii++) {
			if(g1[ii] != 0 && g1[ii-1] == 0) {
				gendist_minmax[0] = g0[ii];
				break;
			}
		}
		cerr << "Finding max\n";
		for(int ii = g0.size() - 1; ii > 0; ii--) {
			if(g1[ii] != g1[ii-1]) {
				gendist_minmax[1] = g0[ii];
			}
		}
		cerr << "Making distances\n";
		for(int i = 0; i < snp_count; i++) {
			int j;
			if(snps[i].position <= g0[0]) {
				double startrate = (g1[1]-g1[0])/(g0[1]-g0[0]);
				snps[i].genetic_position = startrate*(snps[i].position-g0[0]);
			} else {
				for(j = 0; j < g0.size()-1; j++) {
					if(snps[i].position <= g0[j]) {
						snps[i].genetic_position = g1[0];
						break;
					}
					if((snps[i].position >= g0[j]) && (snps[i].position <= g0[j+1])) {
						snps[i].genetic_position = g1[j]+(g1[j+1]-g1[j])/(g0[j+1]-g0[j])*(snps[i].position-g0[j]);
						break;
					}
				}
				if(j == g0.size() - 1) {
					double endrate = (g1[g1.size()-1]-g1[g1.size()-2])/(g0[g0.size()-1]-g0[g0.size()-2]);
					snps[i].genetic_position = snps[i-1].genetic_position+endrate*(snps[i].position - snps[i-1].position);
				}
			}
		}
	}
};

list<pint> createSubPairs(int ac) {
	vector<pint> hold;
	list<pint> out;
	for (int i = 0; i < ac; i++) {
		for (int j = 0; j < i; j++) {
			pint p;
			p.first = i;
			p.second = j;
			hold.push_back(p);
		}
	}
	random_shuffle(hold.begin(),hold.end());
	for (int i = 0; i < 5*ac; i++) {
		out.push_back(hold[i]);
	}
	return out;
}

list<pint> createPairs(int ac) {
	list<pint> out;
	for(int i = 0; i < ac; i++) {
		for (int j = 0; j < i; j++) {
			pint p;
			p.first = i;
			p.second = j;
			out.push_back(p);
		}
	}
	return out;
}


class rare_matrix {
	public:
	int ** distmat_left;
	int ** distmat_right;
	list<pint> active_idx, active_left, active_right;
	int * hap_idx;
	int ac;
	int snp_num;
	rare_matrix() { };
	rare_matrix(int allele_count, vector<int> & rs, vector<int> & rd, int sn) {
		ac = allele_count;
		snp_num = sn;
		cout << "SN: " << sn << " SNP " << snp_num << endl;
		hap_idx = new int[allele_count];
		distmat_left = new int*[allele_count];
		distmat_right = new int*[allele_count];
		if (allele_count > subsample_pairs) {
			active_idx = createSubPairs(allele_count);
		} else {
			active_idx = createPairs(allele_count);
		}
		active_left = active_idx;
		active_right = active_idx;
		for(int i = 0; i < rs.size(); i++) {
			hap_idx[i] = rs[i];
		}
		for(int i = rs.size(); i < rs.size()+rd.size(); i++) {
			hap_idx[i] = rd[i-rs.size()];
		}
		for(int i = 1; i < allele_count; i++) {
			distmat_left[i] = new int[i];
			distmat_right[i] = new int[i];
			for (int j = 0; j < i; j++) {
				distmat_left[i][j] = -2;
				distmat_right[i][j] = -2;
			}
		}
		//for (int i = 0; i < active_idx.size(); i++) {
		for (auto itr = active_idx.begin(); itr != active_idx.end(); itr++) {
			distmat_left[itr->first][itr->second] = -1;
			distmat_right[itr->first][itr->second] = -1;
		}
	}
	bool fillMat(vcf_data & vcf, int snp_idx, int ** distmat, list<pint> & idx) {
		//int total_left = ac*(ac-1)/2;
		/*for (int i = 1; i < ac; i++) {
			for (int j = 0; j < i; j++) {
				if (distmat[i][j] != -1) {
					//total_left -= 1;
					continue;
				}
				char hap1 = vcf.getHap(hap_idx[i],snp_idx);
				char hap2 = vcf.getHap(hap_idx[j],snp_idx);
				if (hap1 != hap2) {
					distmat[i][j] = snp_idx;
				}
			}
		}*/
		auto el = idx.begin();
		while (el != idx.end()) {
			int i = el->first;
			int j = el->second;
			if (distmat[i][j] != -1) { 
				el++;				
				continue;
			}
			char hap1 = vcf.getHap(hap_idx[i],snp_idx);
			char hap2 = vcf.getHap(hap_idx[j],snp_idx);
			if (hap1 != hap2) {
				distmat[i][j] = snp_idx;
				idx.erase(el++);
			}
			else { 
				el++;
			}
		}
		return (idx.size() == 0);
					

	}
	void endSide(int val, int ** distmat) {
		for (int i = 1; i < ac; i++) {
			for (int j = 0; j < i; j++) {
				if (distmat[i][j] == -1) {
					distmat[i][j] = val;
				}
			}
		}
	}
	void fillMats(vcf_data & vcf) {
		int i;
		for(i = snp_num-1; i > 0; i--) {
			if (!vcf.snps[i].isValidForComp()) { continue; }
			if (fillMat(vcf,i,distmat_left,active_left)) { break; }
		}
		if (i <= 0) { endSide(0,distmat_left); }
		for (i = snp_num+1; i < vcf.snp_count - 1; i++) {
			if (!vcf.snps[i].isValidForComp()) { continue; }
			if(fillMat(vcf,i,distmat_right,active_right)) { break; }
		}
		if (i >= vcf.snp_count - 1) {
			endSide(vcf.snp_count-1,distmat_right);
		}
		

	}
};

string ast(int a, bool b) {
	cout << a << endl;
	stringstream s("");
	s << a;
	if(b) {
		s << '*';
	}
	s << '\t';
	return s.str();
}

string ast(double a, bool b) {
	cout << a << endl;
	stringstream s("");
	s << a;
	if(b) {
		s << '*';
	}
	s << '\t';
	return s.str();
}


class output_data {
	public:
	int pos, count;
	int p1_5_bp, p2_5_bp, p1_3_bp, p2_3_bp, p1_total_bp, p2_total_bp;
	double p1_5_gen, p2_5_gen, p1_3_gen, p2_3_gen, p1_total_gen, p2_total_gen;
	bool same_1,same_2;
	bool end_flags[12] = {false};
	int single_phase;
	string ann;
	output_data(vcf_data & vcf, hapset & h1, hapset & h2) {
		int front_end = 0;
		int back_end = vcf.snp_count-1;
		single_phase = 0;
		count = vcf.snps[h1.snp_id].allele_count;
		ann = vcf.snps[h1.snp_id].annotation;
		pos = vcf.snps[h1.snp_id].position;
		p1_5_bp = vcf.getSnpDist(h1.data_right.maxhap,h1.snp_id);
		if(h1.data_right.maxhap == back_end) { end_flags[0] = 1; end_flags[4] = 1; }
		p1_3_bp = vcf.getSnpDist(h1.snp_id,h1.data_left.maxhap);
		if(h1.data_left.maxhap == front_end) { end_flags[1] = 1; end_flags[5] = 1; }
		p2_5_bp = vcf.getSnpDist(h2.data_right.maxhap,h2.snp_id);
		if(h2.data_right.maxhap == back_end) { end_flags[2] = 1; end_flags[6] = 1; }
		p2_3_bp = vcf.getSnpDist(h2.snp_id,h2.data_left.maxhap);
		if(h2.data_left.maxhap == front_end) { end_flags[3] = 1; end_flags[7] = 1; }
		p1_5_gen = vcf.getSnpGm(h1.data_right.maxhap,h1.snp_id);
		p1_3_gen = vcf.getSnpGm(h1.snp_id,h1.data_left.maxhap);
		p2_5_gen = vcf.getSnpGm(h2.data_right.maxhap,h2.snp_id);
		p2_3_gen = vcf.getSnpGm(h2.snp_id,h2.data_left.maxhap);
		p1_total_bp = 0, p2_total_bp = 0, p1_total_gen = 0, p2_total_gen = 0;
	}
	output_data(vcf_data & vcf, hapset & h1) {
		int front_end = 0;
		int back_end = vcf.snp_count-1;
		single_phase = 1;
		count = vcf.snps[h1.snp_id].allele_count;
		ann = vcf.snps[h1.snp_id].annotation;
		pos = vcf.snps[h1.snp_id].position;
		p1_5_bp = vcf.getSnpDist(h1.data_right.maxhap,h1.snp_id);
		if(h1.data_right.maxhap == back_end) { end_flags[0] = 1; end_flags[4] = 1; }
		p1_3_bp = vcf.getSnpDist(h1.snp_id,h1.data_left.maxhap);
		if(h1.data_left.maxhap == front_end) { end_flags[1] = 1; end_flags[5] = 1; }
		p2_5_bp = 0;
		p2_3_bp = 0;
		p1_5_gen = vcf.getSnpGm(h1.data_right.maxhap,h1.snp_id);
		p1_3_gen = vcf.getSnpGm(h1.snp_id,h1.data_left.maxhap);
		p2_5_gen = 0;
		p2_3_gen = 0;

	}
	output_data() { }
	void printData(bool mult) {
		*of_vec[count - range_start] << chr_num << '\t' << pos << '\t';
		*of_vec[count - range_start] << ast(p1_5_bp,end_flags[0]);
		*of_vec[count-range_start] << ast(p1_3_bp,end_flags[1]) << ast(p2_5_bp,end_flags[2]) << ast(p2_3_bp,end_flags[3]);
		if(no_gen_map == 0) { *of_vec[count-range_start] << ast(p1_5_gen,end_flags[4]) << ast(p1_3_gen,end_flags[5]) << ast(p2_5_gen,end_flags[6]) << ast(p2_3_gen,end_flags[7]); }
		*of_vec[count-range_start] << ann;
		if(!mult) { *of_vec[count-range_start] << endl; }
	}



};



class mx_output {
	public:
	list<output_data> od;
	rare_matrix mat;
	int ac;
	int snp_idx;
	mx_output() { };
	mx_output(list<output_data> & l, rare_matrix & m) {
		od = l;
		mat = m;
		ac = m.ac;
		snp_idx = m.snp_num;
		cout << "M: " << m.snp_num << " SI " << snp_idx << endl;
	}
	void print() {
		cout << "AC: " << ac << endl;
		cout << "SNP_IDX: " << snp_idx << endl;
		for(int i = 0; i < ac; i++) {
			for(int j = 0; j < ac; j++) {
				if (i > j) { 
					cout << mat.distmat_left[i][j] << '\t';
				} else if (i < j) {
					cout << mat.distmat_left[j][i] << '\t';
				} else { cout << "0\t"; }
			} 
			cout << "\t";
			for(int j = 0; j < ac; j++) {
				if (i > j) { 
					cout << mat.distmat_right[i][j] << '\t';
				} else if (i < j) {
					cout << mat.distmat_right[j][i] << '\t';
				} else { cout << "0\t"; }
			}
			cout << endl;
		}
	}
		
};

int ** intInit(int s) {
	int ** out = new int*[s];
	for(int i = 0; i < s; i++) {
		out[i] = new int[s];
		for(int j = 0; j < s; j++) {
			out[i][j] = 0;
		}
	}
	return out;
}

double ** doubleInit(int s) {
	double ** out = new double*[s];
	for(int i = 0; i < s; i++) {
		out[i] = new double[s];
		for(int j = 0; j < s; j++) {
			out[i][j] = 0;
		}
	}
	return out;
}

void boolInit(bool ** b, int s) {
	cout << "BI " << s << endl;
	b = new bool*[s];
	for (int i = 0; i < s; i++) {
		b[i] = new bool[s];
		for(int j = 0; j < s; j++) {
			b[i][j] = 0;
		}
	}
}

class mx_print{
	public:
	int ** left_pos;
	int ** right_pos;
	double ** left_gen;
	double ** right_gen;
	int ** left_ends;
	int ** right_ends;
	int ac;
	int snp_idx;
	mx_print() { };
	mx_print(vcf_data & vcf, mx_output & m) {
		ac = m.ac;
		snp_idx = m.snp_idx;
		left_pos = intInit(ac);
		right_pos = intInit(ac);
		left_gen = doubleInit(ac);
		right_gen = doubleInit(ac);
		left_ends = intInit(ac);
		right_ends = intInit(ac);
		for (int i = 0; i < ac; i++) {
			for (int j = 0; j < ac; j++) {
				if (i == j) {
					left_pos[i][j] = 0;
					right_pos[i][j] = 0;
					left_gen[i][j] = 0;
					right_gen[i][j] = 0;
					left_ends[i][j] = 0;
					right_ends[i][j] = 0;
				} else if (i > j) {
					right_pos[i][j] = vcf.getSnpDist(m.mat.distmat_right[i][j],snp_idx);
					left_pos[i][j] = vcf.getSnpDist(snp_idx,m.mat.distmat_left[i][j]);
					right_gen[i][j] = vcf.getSnpGm(m.mat.distmat_right[i][j],snp_idx);
					left_gen[i][j] = vcf.getSnpGm(snp_idx,m.mat.distmat_left[i][j]);
					if (m.mat.distmat_right[i][j] == vcf.snp_count-1) { right_ends[i][j] = 1; }
					if (m.mat.distmat_left[i][j] == 0) { left_ends[i][j] = 1; }
				} else {
					right_pos[i][j] = vcf.getSnpDist(m.mat.distmat_right[j][i],snp_idx);
					left_pos[i][j] = vcf.getSnpDist(snp_idx,m.mat.distmat_left[j][i]);
					right_gen[i][j] = vcf.getSnpGm(m.mat.distmat_right[j][i],snp_idx);
					left_gen[i][j] = vcf.getSnpGm(snp_idx,m.mat.distmat_left[j][i]);
					if (m.mat.distmat_right[j][i] == vcf.snp_count-1) { 
						right_ends[i][j] = 1; }
					if (m.mat.distmat_left[j][i] == 0) { 
						left_ends[i][j] = 1; }

				}
			}
		}
	}
	void printRow(int row) {
		*of_vec[ac-range_start] << '\t';
		for(int j = 0; j < ac; j++) {
			*of_vec[ac-range_start] << ast(right_pos[row][j],right_ends[row][j]);
		}
		for(int j = 0; j < ac; j++) {
			*of_vec[ac-range_start] << ast(left_pos[row][j],left_ends[row][j]);
		}
		if(no_gen_map == 0) {
			for(int j = 0; j < ac; j++) {
				*of_vec[ac-range_start] << ast(right_gen[row][j],right_ends[row][j]);
			}
			for(int j = 0; j < ac; j++) {
				*of_vec[ac-range_start] << ast(left_gen[row][j],left_ends[row][j]);
			}
		}
		*of_vec[ac-range_start] << endl;
	}
		

};

uchar * compressHapRow(vector<char> & v, int snp_count) {
	int comp_size = (snp_count+7)/8;
	//cout << comp_size << endl;
	uchar * out = new uchar[comp_size];
	for(int i = 0; i < comp_size; i++) {
		int t = 0;
		for(int j = 0; j < 8; j++) {
			int reg = pow(2,7-j);
			int old_idx = 8*i+j;
			if(old_idx >= snp_count) { break; }
			t += reg*int(v[old_idx]);
			//cout << t << " " << int(v[old_idx]) << endl;
		}
		out[i] = (unsigned char)t;
		//cout << int(out[i]) << endl;
	}
	return out;
}


uchar ** convertVecToData(vvchar & hap_hold) {
	cout << compress_flag << endl;
	int hap_count = hap_hold.size();
	int snp_count = hap_hold[0].size();
	uchar ** d = new uchar*[hap_count];
	for(int i = 0; i < hap_hold.size(); i++) {
		/*if(compress_flag == 1) {
			d[i] = compressHapRow(hap_hold[i],snp_count);
		} else {*/
			//cout << snp_count << endl;
		d[i] = new uchar[snp_count];
		for(int j = 0; j < snp_count; j++) {
			d[i][j] = hap_hold[i][j];
		}
		//}
		//hap_hold[i].clear();
		(vector<char>()).swap(hap_hold[i]);
	}
	(vvchar()).swap(hap_hold);
	cout << "After swaps\n";
	return d;


}




vector<string> pullNames(string line) {
        stringstream s(line);
        string hold, junk;
        s >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk; //9 junks, verify that merged files didn't skip a person
        vector<string> names;
        while(s >> hold) {
                names.push_back(hold);
        }
        return names;
}

double getGenCol(string line, int col, double & phys) {
	stringstream s(line);
	s >> phys;
	double j = 0, out = 0;
	for(int i = 0; i < col; i++) {
		s >> j;
	}
	bool b = (s >> out);
	if(!b) { return -1.0; }
	return out;
}

list<int> subsampleIdx(vector<int> & vals, int subsample) {
	random_shuffle(vals.begin(), vals.end());
	list<int> out;
	for(int i = 0; i < subsample; i++) {
		out.push_back(vals[i]);
	}
	out.sort();
	return out;
}


list<int> getSampleListSingle(vcf_data & vcf, int idx, vector<int> & rare_idx, int subsample) {
	list<int> out; 
	vector<int> vhold;
	for(int i = 0; i < vcf.hap_count; i+=2) {
		char h1 = vcf.getHap(i,idx);
		char h2 = vcf.getHap(i+1,idx);
		if(h1 == 0 && h2 == 0) {
			vhold.push_back(i);
			vhold.push_back(i+1);
		} else if (h1 == 1) {
			rare_idx.push_back(i);
			rare_idx.push_back(i+1);
		} else {
			rare_idx.push_back(i+1);
			rare_idx.push_back(i);
		}
	}
	if(subsample == 0) {
		//cout << "Adding to list\n";
		for(int i = 0; i < vhold.size(); i++) {
			out.push_back(vhold[i]);
		}
	} else {
		out = subsampleIdx(vhold,subsample);
	}
	return out;
}

//For subsampling, 0 = dont, anything positive means list will contain subsample number of indices
list<int> getSampleListMult(vcf_data & vcf, int idx, vector<int> & rare_single, vector<int> & rare_double, int subsample) {
	list<int> out;
	vector<int> vhold;
	for(int i = 0; i < vcf.hap_count; i+=2) {
		char h1 = vcf.getHap(i,idx);
		char h2 = vcf.getHap(i+1,idx);
		if(h1 == 0 && h2 == 0) {
			vhold.push_back(i);
			vhold.push_back(i+1);
		} else if(h1 == 1 && h2 == 1) {
			rare_double.push_back(i);
			rare_double.push_back(i+1);
		} else if(h1 == 1) {
			rare_single.push_back(i);
			vhold.push_back(i+1);
		} else {
			rare_single.push_back(i+1);
			vhold.push_back(i);
		}
		
	}
	//Subsample here
	if(subsample == 0) {
		//cout << "Adding to list\n";
		for(int i = 0; i < vhold.size(); i++) {
			out.push_back(vhold[i]);
		}
	} else {
		out = subsampleIdx(vhold,subsample);
	}

	return out;
}

list<int> getModList(list<int> & haps, int sw) {
	list<int> out;
	for(auto i = haps.begin(); i != haps.end(); ++i) {
		if(*i != sw) { out.push_back(*i); }
		else {
			if(sw%2 == 1) {
				out.push_back(*i-1);
			} else {
				out.push_back(*i+1);
			}
			
		}
	}
	return out;

}

bool fillStoppingPoint(vcf_data & vcf, hapset & h, maxhapdata & mhd, list<int> & comp_haps, int snp_idx) {
	list<int>::iterator i = comp_haps.begin();
	while(i != comp_haps.end()) {
		if(mhd.cutoffs[*i] != -1) {
			comp_haps.erase(i++);
			continue;
		}
		char rare_hap = vcf.getHap(h.hap_id,snp_idx);
		char current_hap = vcf.getHap(*i,snp_idx);
		if(rare_hap != current_hap) {
			mhd.cutoffs[*i] = snp_idx;
			comp_haps.erase(i++);
			mhd.remaining -= 1;
		} else { 
			++i;
		}

	}
	if(comp_haps.size() == 0) { return 1; } //Return 1 means loop in calling function should end
	return 0;
}

bool fillStoppingRare(vcf_data & vcf, hapset & h, maxhapdata & mhd, int snp_idx) {
	bool included = 0;
	for(auto ii = vcf.snps[snp_idx].hap_ids.begin(); ii != vcf.snps[snp_idx].hap_ids.end(); ii++) {
		//cout << *ii << '\t';
		if (h.hap_id == *ii) {
			included = 1;
			break;
		}
	}
	//cout << endl;
	if (included) { return 0; }
	for(auto ii = vcf.snps[snp_idx].hap_ids.begin(); ii != vcf.snps[snp_idx].hap_ids.end(); ii++) {
		//cout << *ii << "\t" << mhd.cutoffs[*ii] << "\n";
		if(mhd.cutoffs[*ii] == -1) {
			mhd.cutoffs[*ii] = snp_idx;
			mhd.remaining -= 1;
			//cout << "found: " << mhd.cutoffs[*ii] << endl;
		}
	}
	return 1;
}


void findMax(vcf_data & vcf, hapset & h) {
	list<int> c_left = h.hap_idx_list;
	//left portion
	int i = h.snp_id-1;
	for(i = h.snp_id-1; i > 0; i--) {
		if(!vcf.snps[i].isValidForComp()) { continue; }
		if(vcf.snps[i].rare) {
			bool success = fillStoppingRare(vcf,h,h.data_left,i);
			if (h.data_left.remaining == 0) { 
				//cout << "ENDED\n";
				break; }
			if (success) { continue; }
		}
		bool filled = fillStoppingPoint(vcf,h,h.data_left,c_left,i);
		if(filled) {
			break;
		}
	}
	h.data_left.maxhap = (i >= 0 ? i : 0);
	//cout << i << endl;
	if(i <= 0) {
		h.fillLeft(0);
		//cout << "filling left side\n";
	}
	list<int> c_right = h.hap_idx_list;
	i = h.snp_id+1;
	for(i = h.snp_id+1; i < vcf.snp_count - 1; i++) {
		if(!vcf.snps[i].isValidForComp()) { continue; }
		if(vcf.snps[i].rare) {
			bool success = fillStoppingRare(vcf,h,h.data_right,i);
			if (h.data_right.remaining == 0) { 
				//cout << "ENDED\n";
				break; }
			if (success) { continue; }
		}
		bool filled = fillStoppingPoint(vcf,h,h.data_right,c_right,i);
		if(filled) {
			break;
		}
	}
	h.data_right.maxhap = i;
	//cout << i << " " << vcf.snp_count << endl;
	if(i >= vcf.snp_count-1) {
		h.fillRight(vcf.snp_count-1);
		h.data_right.maxhap = vcf.snp_count-1;
		//cout << "filling right side " << vcf.snp_count << endl;
	}

}

void printVerify(vcf_data & vcf, hapset & h) {
	cout << "HAP ID: " << h.snp_id << endl;
	cout << h.hap_id << '\t';
	for(int i = h.snp_id; i >= h.data_left.maxhap; i--) {
		cout << int(vcf.getHap(h.hap_id,i));
	}
	cout << endl;
	for(auto itr = h.hap_idx_list.begin(); itr != h.hap_idx_list.end(); ++itr) {
		cout << *itr << '\t' << h.data_left.cutoffs[*itr] << '\t';
		for(int i = h.snp_id; i >= h.data_left.cutoffs[*itr]; i--) {
			cout << int(vcf.getHap(*itr,i));
		}
		cout << endl;
	}
	cout << endl;
	cout << h.hap_id << '\t';
	for(int i = h.snp_id; i <= h.data_right.maxhap; i++) {
		cout << int(vcf.getHap(h.hap_id,i));
	}
	cout << endl;
	for(auto itr = h.hap_idx_list.begin(); itr != h.hap_idx_list.end(); ++itr) {
		cout << *itr << '\t' << h.data_right.cutoffs[*itr] << '\t';
		for(int i = h.snp_id; i <= h.data_right.cutoffs[*itr]; i++) {
			cout << int(vcf.getHap(*itr,i));
		}
		cout << endl;
	}
}


mx_output findLongestHaps(vcf_data & vcf, int idx) {
	mx_output out;
	vector<int> rare_idx;
	vector<hapset> vh;
	list<output_data> out_list;
	out.ac = vcf.snps[idx].allele_count;
	if (matrix_flag != 1 && vcf.snps[idx].allele_count == 1) {
	//if(vcf.snps[idx].allele_count == 1) {
		vector<int> rare_idx;
		list<int> comparison_haps = getSampleListSingle(vcf,idx,rare_idx,0);
		hapset h1(rare_idx[0],idx,vcf.hap_count,comparison_haps,rare_idx[1]);
		hapset h2(rare_idx[1],idx,vcf.hap_count,comparison_haps,rare_idx[0]);
		findMax(vcf,h1);
		findMax(vcf,h2);
		output_data od(vcf,h1,h2);
		out.od.push_back(od);
	} else { 
		vector<int> rare_single, rare_double;
		list<int> comparison_haps = getSampleListMult(vcf,idx,rare_single,rare_double,0);
		if (matrix_flag != 1) {
		cout << "GENERATING STUFF\n";
		for(int i = 0; i < rare_single.size(); i++) {
			hapset h1(rare_single[i],idx,vcf.hap_count,comparison_haps);
			int sw_person;
			if(rare_single[i]%2==1) {
				sw_person = rare_single[i]-1;
			} else {
				sw_person = rare_single[i]+1;
			}
			list<int> swap_phase_haps = getModList(comparison_haps,sw_person);
			hapset h2(sw_person,idx,vcf.hap_count,swap_phase_haps);
			findMax(vcf,h1);
			findMax(vcf,h2);
			output_data od(vcf,h1,h2);
			out.od.push_back(od);
			//od.printData();
		}
		for(int i = 0; i < rare_double.size(); i+=2) {
			hapset h1(rare_double[i],idx,vcf.hap_count,comparison_haps);
			findMax(vcf,h1);
			output_data od(vcf,h1);
			out.od.push_back(od);
			//od.printData();
			hapset h2(rare_double[i+1],idx,vcf.hap_count,comparison_haps);
			findMax(vcf,h2);
			output_data od2(vcf,h2);
			out.od.push_back(od2);
			//od2.printData();
		}
		}
		if(matrix_flag != 0) {
			rare_matrix rm(vcf.snps[idx].allele_count,rare_single,rare_double,idx);
			rm.fillMats(vcf);
			out.mat = rm;
			out.ac = rm.ac;
			out.snp_idx = rm.snp_num;
			out.print();
		}
	}
	//for (auto ii = out_list.begin(); ii != out_list.end(); ii++) {
	//	out.od.push_back(*ii);
	//}
		
	//out.od = out_list;
	return out;
}
		
int printCurrentResults(vcf_data & vcf, mx_output * ol, bool * flags, int max, int & prevcount) {
	for (int i = prevcount; i < max; i++) {
		if (flags[i] == 0) { return 0; }
		bool mult = (ol[i].ac != 1) && (matrix_flag != 0);
		int cur_idx = 0;
		mx_print *mp;
		if(mult) {
			//mx_print t(vcf,ol[i]);
			mp = new mx_print(vcf,ol[i]);
		}
		for (auto itr = ol[i].od.begin(); itr != ol[i].od.end(); ++itr) {
			if(matrix_flag != 1) {
				(*itr).printData(mult);
			}
			if(mult) {
				mp->printRow(cur_idx);
			}
			cur_idx++;
		}
		prevcount++;
	}
	return 1;
}



int main(int argc, char ** argv) {
	setbuf(stdout,NULL);
	int seed = 1;
	string map_name, vcf_name = "", output_tag = "rare_allele_run";
	vcf_data v;
	int range_end = 10;
	int pos_start = 0;
	for(int i = 1; i < argc; i++) {
		string arg = argv[i];
		if(arg == "-i") { vcf_name = argv[++i]; }
		else if(arg == "-m") { map_name = argv[++i]; no_gen_map = 0; }
		else if(arg == "-r") { range_start = atoi(argv[++i]); range_end = atoi(argv[++i]); }
		else if(arg == "-o") { output_tag = argv[++i]; }
		else if(arg == "-c") { compress_flag = 1; }
		else if(arg == "-p") { pos_start = atoi(argv[++i]); }
		else if(arg == "--ignore-singletons") { ignore_singletons = 1; }
		else if(arg == "-matf") { matrix_flag = atoi(argv[++i]); }
		else {
			cerr << "Invalid option: " << arg << endl;
			return 0;
		}

	}
	//cout << ignore_singletons << endl;
	srand(seed);
	output_tag += ".results";
	v.readVcf(vcf_name);
	if (no_gen_map == 0) {
		v.addGenData(map_name,0);
	}
	for (int i = 0; i <= range_end - range_start; i++) {
		string filename = output_tag+"."+to_string(range_start+i);
		ofstream *op = new ofstream(filename.c_str());
		of_vec.push_back(op);
	}
	vector<int> rareIdx = v.getRareIdx(range_start,range_end,pos_start);

	//v.print();
	//list<output_data> * output_hold = new list<output_data>[rareIdx.size()];
	mx_output * output_hold = new mx_output[rareIdx.size()];
	bool * output_flag = new bool[rareIdx.size()];

	for(int i = 0; i < rareIdx.size(); i++) {
		output_flag[i] = 0;

	}
	//return 0;
	int output_loc = 0;
	#pragma omp parallel for ordered schedule(dynamic,50)
	for (int ii = 0; ii < rareIdx.size(); ii++) {
		int i = rareIdx[ii];
		output_hold[ii] = findLongestHaps(v,i);
		output_flag[ii] = 1;
		if (ii%10 == 0) {
			#pragma omp critical
			printCurrentResults(v,output_hold,output_flag,rareIdx.size(),output_loc);
		}
	}
	printCurrentResults(v,output_hold,output_flag,rareIdx.size(),output_loc);
	return 0;
}


















