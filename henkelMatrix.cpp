#include <iostream>
#include <string>
#include <math.h>
#include <cmath>
#include <boost/multi_array.hpp>

//the code first compute the left superior triangle (i<m, j<n-m) and then the
//second inferior triangle (i<m j>n-m). The coefficients are computed on the first
//an last line of the matrix, and then we move along the diagonal to fill the matrix

//enventually, the matrix is printed per line

void transfo_henkel(boost::multi_array<int,2>  &input){
	int m =int(input.shape()[0]);
	int n=int(input.shape()[1]);
	if(m>n){std::cout << "Error from the data:" << std::endl;}
	else{
		boost::multi_array<float,2> output(boost::extents[m][n]);
		for(int j=0;j<n;j++){
			int l= j;
			int k=0;
			int my_cont=0;
			for(int i=0;i<m;i++){if (l>=0){k=k+input[i][l];l=l-1;my_cont++;}
								 else{break;}}
			output[0][j]=float(k/my_cont);

		}
		for(int j=n-m;j<n;j++){
			int l=j;
			int k=0;
			int my_cont=0;
			for(int i=m-1;i>=0;i--){if(l<n){k=k+input[i][l];l=l+1;my_cont++;}
									else{break;}}
				output[m-1][j]=float(k/my_cont);


		}
		
	for(int j=1;j<n;j++){
		int l=j-1;
		for(int i=1;i<m;i++){if (l>=0){output[i][l]=output[0][j];l=l-1;}
								 else{break;}}


	}

	for(int j=n-m;j<n;j++){
		int l=j+1;
		for(int i=m-2;i>=0;i--){if(l<=n-1){output[i][l]=output[m-1][j];l=l+1;}
								else{break;}}
	}
	
	for(int i=0;i<m;i++){
		std::cout<<"line "<<i<<std::endl;
		for(int j=0;j<n;j++){
			std::cout << output[i][j] << std::endl;
		}
	}

}
}
