

//#redefinde stdout printing for better user feedback

using namespace std;

#define lL " \E[37;45m\033[1m>\033[0m "
#define wL " \E[37;31m\033[1m>\033[0m "
#define cL " \E[37;31m\033[1m   \033[0m "
#define printL cout<<lL<<
#define printW cout<<wL<<
#define printC cout<<cL<<
#define endL <<endl
#define emptyLine (cout<<endl)
#define getStr() .str().c_str()
