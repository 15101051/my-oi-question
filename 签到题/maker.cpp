#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<ctime>
#include<cctype>
#include<cstring>
using namespace std;
inline int rd()
{
	register int res = 0,f = 1;register char c = getchar();
	while(!isdigit(c)){if(c == '-')f = -1;c = getchar();}
	while(isdigit(c))res = (res << 1) + (res << 3) + c - '0',c = getchar();
	return res * f;
}
int main()
{
	freopen("sign9.in","w",stdout);
	srand(time(NULL));
	cout << 10000000000 << " " << 1ll * rand() * rand() % 1000000000 + 1 << endl;
	
	return 0;
}

