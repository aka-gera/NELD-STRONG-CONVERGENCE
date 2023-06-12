//-------------------------------------- 
//   Functions required in 'input.hpp' 
//--------------------------------------

#include"input.hpp"
//--------- start reading of the input file -------------
void input::load(void){

//------------- Reading the file line by line --------------------
	ifstream from("param.in");
	read_item<double>(from,"Time step                 (dt) ",&t_step);
	read_item<double>(from,"Thermalization time (T)        ",&T_therm);
	read_item<double>(from,"Simulation time (T)            ",&T_final);
	read_item<int>   (from,"Number of Xmakemol outputs     ",&n_xmakemol);
	read_item<int>   (from,"Number of Observable outputs   ",&n_obs);
	read_item<double>(from,"Friction coefficient      (xi) ",&xi);
	read_item<double>(from,"Inverse temperature     (beta) ",&beta);
	read_item<int>   (from,"Particles in x        (Ndim_x) ",&Ndim_x);
	read_item<int>   (from,"Particles in y        (Ndim_y) ",&Ndim_y);
	read_item<int>   (from,"Particles in z        (Ndim_z) ",&Ndim_z);
	read_item<double>(from,"Elementary cell size       (a) ",&a);
	read_item<double>(from,"WCA equilibrium distance (sig) ",&sig);
	read_item<double>(from,"WCA energy epsilon       (eps) ",&eps);
	read_item<double>(from,"Cut off radius   (sigma units) ",&dcut);
	read_item<Matrix>(from,"Matrix A  (Diagonalized)       ",&A);
	read_item<Matrix>(from,"Matrix S                       ",&S);
	read_item<Matrix>(from,"Matrix S^-1                    ",&Sinv);
	read_item<bool>  (from,"Use GenKR                      ",&genKR);

	dim = 3;
}

template <class T>

//---------- Reading a line in the input file -----------------
void input::read_item(istream& from,const char *_title, T * val){
	char title[80];
	char c;

	from.get(title,80,':');
	std::streamsize size_read = from.gcount();
	from.get(c);
	if(strcmp(title,_title)){ 
		cerr << "The item reads '" << title << "' different from expected.'"
			<< endl;
		cout << "Input " << _title << endl;
		cin >> *val; 
		from.seekg(-(size_read+1), ios::cur);
	}else{
		from >> *val;
		from.get(title,80,'\n');
		from.get(c);
	}
}





