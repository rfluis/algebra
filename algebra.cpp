#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <iomanip>

namespace rcl{
  int strlen(const char* str){int idx=0; while(str[idx++]); return idx;}
  void memcpy(char *dest,const char *src,int bytes){for (int idx=0;idx<bytes;idx++) dest[idx]=src[idx]; return;} 
  class String{
    char *ptr=NULL;
    int sz=0;
  public:
    void fit(){ptr=(char*)realloc(ptr,sz); return;}
    void fit(const int& n){sz=n; fit(); return;}
    String(){return;}
    String(const char* str){fit(rcl::strlen(str));rcl::memcpy(ptr,str,sz); return;}
    String(const int& n,const char& c){fit(n); for (int idx=0;idx<n;idx++) ptr[idx]=c; return;}
    int size()const{return sz;}
    void push_back(const char& c){fit(sz+1);ptr[sz-1]=c; return;}
    char& operator [](const int& idx){return ptr[idx];}
    const char& operator [](const int& idx)const{return ptr[idx];}
    };
  }

void hbar(){std::cout << std::string(80,'=') << std::endl;}
std::string helper(int width, const std::string& str) {
    int len = str.length();
    if(width < len) { return str; }

    int diff = width - len;
    int pad1 = diff/2;
    int pad2 = diff - pad1;
    return std::string(pad1, ' ') + str + std::string(pad2, ' ');
}

std::string translate(const std::string& src,const char& key,const char& patch){
  std::string ans;
  for (std::string::const_iterator itr=src.begin();itr!=src.end();itr++){
    if (*itr==key) ans.push_back(patch);
    else ans.push_back(*itr);
    }
  return ans;
  }

std::string producto(const std::string& a,const std::string& b){
  if (a=="0") return "0";
  if (b=="0") return "0";
  if (a=="1") return b;
  if (b=="1") return a;
  //std::cout << "Producto:" << a << " * " << b << std::endl;
  return (a+"*"+b);
  }

class expr_product;
std::vector<expr_product> parse_eqn(const std::string& eq);

template<typename T> class matrix{
  int xx,yy;
  std::vector<T> data;
  public:
  matrix(){xx=yy=0; return;}
  matrix(const int& x,const int& y){fit(x,y);return;}
  void fit(){data.resize(xx*yy); return;}
  void fit(const int& x,const int& y){xx=x;yy=y;fit();return;}
  int gx()const{return xx;}
  int gy()const{return yy;}
  std::string log(){
    std::stringstream ss;
    ss << "(" << gx() << "," << gy() << ":" << data.size() << ")";
    return ss.str();
    }
  T& operator()(const int& x,const int& y){return data[x+y*xx];}
  const T& operator()(const int& x,const int& y)const {return data[x+y*xx];}
  T& operator[](const int& x){return data[x];}
  const T& operator[](const int& x)const{return data[x];}
  matrix transpose(){
    matrix ans(yy,xx);
    for (unsigned int idy=0;idy<yy;idy++) for (unsigned int idx=0;idx<xx;idx++) ans(idy,idx)=(*this)(idx,idy);
    return ans;
    }
  matrix reverse(){
    matrix ans(xx,yy);
    for (unsigned int idy=0;idy<yy;idy++) for (unsigned int idx=0;idx<xx;idx++) ans(idx,idy)=(*this)(xx-1-idx,yy-1-idy);
    return ans;
    }
  std::string LaTeX()const{
    std::stringstream ss;
    ss << "\\left( " << std::endl;
    ss << "\\begin{array}{"<< std::string(xx,'c') <<"}" << std::endl;
    for (int idy=0;idy<yy;idy++){
      for (int idx=0;idx<xx;idx++){
        if (idx!=0) ss << "& ";
        ss << data[idx+idy*xx] << " ";
        }
      if (idy!=yy-1) ss << "\\\\";
      ss << std::endl;
      }
    ss << " \\end{array}" << std::endl;
    ss << "\\right)" << std::endl;
    return ss.str();
    }  
};

matrix<std::string> operator*(const matrix<std::string>& lhs,const matrix<std::string>& rhs){
  if (lhs.gx()!=rhs.gy()) return matrix<std::string>(1,1);
  int nx=rhs.gx(),ny=lhs.gy(),nc=lhs.gx();
  matrix<std::string> p(nx,ny);
  for (unsigned int idy=0;idy<ny;idy++)
  for (unsigned int idx=0;idx<nx;idx++){
    std::string fcproducto;
    for (unsigned int idz=0;idz<nc;idz++){
      std::string subproducto = producto(lhs(idz,idy),rhs(idx,idz));
      if (subproducto!="0") {
        if (!fcproducto.empty()) fcproducto+=" + ";
        fcproducto += subproducto;
      }
    }
   p(idx,idy)=fcproducto;
 }
 return p;
}


void build_lu_matrix_indeces(std::stringstream& ss,const unsigned int& idx,const unsigned int& idy,const bool& latex){
  if (latex) ss << "_{{";
  ss << idy+1 ;
  if (latex) ss << "}{";
  ss << idx+1 ;
  if (latex) ss << "}}";
  return;
  }
void build_lu_matrix_indeces(std::stringstream& ss,const unsigned int& idx,const bool& latex){
  if (latex) ss << "_{";
  ss << idx+1 ;
  if (latex) ss << "}";
  return;
  }

matrix<std::string> build_row_matrix(const int& n,const std::string& kindfmt){
  matrix<std::string> ans(1,n);
  std::string kindstr;
  std::string kindtex;
  bool latex_output=false;
  std::stringstream kindss(kindfmt);
  kindss >> kindstr;
  if (kindss.good()) kindss >> kindtex;
  if (kindtex!="") latex_output=true;
  
  for (unsigned int idx=0;idx<n;idx++){
    std::stringstream ss;
    ss << kindstr;
    build_lu_matrix_indeces(ss,idx,latex_output);
    ans(0,idx)=ss.str();
    }
  return ans;
  }
  
matrix<std::string> build_lu_matrix(const int& n,const std::string& kindfmt){
  enum kind_of_matrix{UNIT_U,UNIT_L,DET_U,DET_L,COMPACT_UNIT_U,COMPACT_UNIT_L,A_MATRIX};
  kind_of_matrix kind;
  std::string kindstr;
  std::string kindtex;
  bool latex_output=false;
  std::stringstream kindss(kindfmt);
  kindss >> kindstr;
  if (kindss.good()) kindss >> kindtex;
  if (kindtex!="") latex_output=true;
  kind=A_MATRIX;
  if (kindstr=="u") kind=UNIT_U;
  if (kindstr=="l") kind=UNIT_L;
  if (kindstr=="U") kind=DET_U;
  if (kindstr=="L") kind=DET_L;
  if (kindstr=="Lu") kind=COMPACT_UNIT_U;
  if (kindstr=="lU") kind=COMPACT_UNIT_L;
  if (kindstr=="uL") kind=COMPACT_UNIT_U;
  if (kindstr=="Ul") kind=COMPACT_UNIT_L;
  matrix<std::string> ans(n,n);
  
  for (unsigned int idy=0;idy<n;idy++)
  for (unsigned int idx=0;idx<n;idx++){
    std::stringstream ss;
    switch(kind){
      case UNIT_U:
        if (idx>idy) ss << "u";
        break;
      case DET_U:
        if (idx>=idy) ss << "u";
        break;
      case UNIT_L:
        if (idx<idy) ss << "l";
        break;
      case DET_L:
        if (idx<=idy) ss << "l";
        break;
      case COMPACT_UNIT_U:
        if (idx<=idy) ss << "l";
        if (idx>idy) ss << "u";
        break;
      case COMPACT_UNIT_L:
        if (idx<idy) ss << "l";
        if (idx>=idy) ss << "u";
        break;
      case A_MATRIX:
      default:
        ss << kindstr;
        break;        
      };
    switch(kind){
      case UNIT_U:
        if (idx==idy) ss << "1";
        else{
          if (idx>idy) build_lu_matrix_indeces(ss,idx,idy,latex_output);
          else ss << "0";
          }
        break;
      case DET_U:
        if (idx>=idy) build_lu_matrix_indeces(ss,idx,idy,latex_output);
        else ss << "0";
        break;
      case UNIT_L:
        if (idx==idy) ss << "1";
        else{
          if (idx<idy) build_lu_matrix_indeces(ss,idx,idy,latex_output);
          else ss << "0";
        }
        break;
      case DET_L:
        if (idx<=idy) build_lu_matrix_indeces(ss,idx,idy,latex_output);
        else ss << "0";
        break;
      case COMPACT_UNIT_U:
        build_lu_matrix_indeces(ss,idx,idy,latex_output);
        break;
      case COMPACT_UNIT_L:
        build_lu_matrix_indeces(ss,idx,idy,latex_output);
        break;
      case A_MATRIX:
      default:
        build_lu_matrix_indeces(ss,idx,idy,latex_output);
        break;        
      };
    ans(idx,idy)=ss.str();
    }
  return ans;
}

class expr_product{
  std::string a,b;
  bool sign; // 0=+ , 1=-
public:
  expr_product(){sign=0;return;}
  expr_product(const bool& s,const std::string& aa,const std::string& bb){a=aa;b=bb,sign=s;return;}
  expr_product operator-(){return expr_product(!sign,a,b);}
  expr_product(const std::string& str){parse_str(str);return;}
  bool is_zero()const{
    if (a=="0"||b=="0"||a==" "||b==" "||a.empty()||b.empty()) return true;
    return false;
    }
  bool is_one()const{
    if (a=="1"&&b=="1"&&sign==0) return true;
    return false;
    }
  bool is_negative()const{return sign;}
  std::string& GetA(){return a;}
  const std::string& GetA()const{return a;}
  std::string& GetB(){return b;}
  const std::string& GetB()const{return b;}
  void exchange(){std::string x=a; a=b;b=x; return;}
  void flip(){sign=!sign; return;}
  void parse_str(const std::string& str){
    sign=false;
    if (str==std::string()) {a="0";b="0";sign=0;}
    std::string::const_iterator prod_operator;
    std::string::const_iterator cursor;
    std::string first,second;
    cursor = std::find(str.begin(),str.end(),'-');
    if (cursor!=str.end()) sign = true;
    prod_operator = std::find(str.begin(),str.end(),'*');
    if (prod_operator==str.end()) second = "1";
    if (cursor==str.end()) cursor=str.begin();
    while ((*cursor==' '||*cursor=='-'||*cursor=='+')&&cursor!=str.end()) cursor++;
    while (*cursor!=' '&&*cursor!='*'&&cursor!=str.end()) first.push_back(*cursor++);
    if (prod_operator!=str.end()){
      cursor = prod_operator;
      while ((*cursor==' '||*cursor=='*')&&cursor!=str.end()) cursor++;
      while (*cursor!=' '&&cursor!=str.end()) second.push_back(*cursor++);
      }
    a=first;
    b=second;
    return;
    }
  std::string str()const{
    std::stringstream ss;
    if (sign)  ss << "- "; else  ss << "+ ";
    ss << a ;
    if (b!="1") ss << " * " << b;
    return ss.str();
    }
  std::string LaTeX()const{
    std::stringstream ss;
    if (is_zero()) return "0";
    if (sign)  ss << "- ";
    ss << a ;
    if (b!="1") ss << " " << b;
    return ss.str();
    }
    std::string LaTeX_sum()const{
    std::stringstream ss;
    if (is_zero()) return "0";
    if (sign)  ss << "- ";
    if (!sign)  ss << "+ ";
    ss << a ;
    if (b!="1") ss << " " << b;
    return ss.str();
    }
  };

std::ostream& operator << (std::ostream& os,const expr_product& e){
  os << e.str();
  return os;
  }

class expr_fraction{
  std::vector<expr_product> num;
  expr_product den;
  public:
  expr_fraction(){return;}
  expr_fraction(const std::string& num,const std::string& den){this->num=parse_eqn(num); this->den=expr_product(den); return;}
  expr_fraction(const std::vector<expr_product>& num,const expr_product& den){this->num=num; this->den=den;return;}
  std::string str()const{
    std::stringstream ss;
    ss << "(";
    for (int idx=0;idx<num.size();idx++) ss << " (" << num[idx] << ") ";
    ss << ") / ( "<< den << " )";
    return ss.str();
    }  
  std::string LaTeX()const{
    std::stringstream ss;
    if (den.is_one()){
      for (unsigned int idx=0;idx<num.size();idx++) {
        if (idx==0) ss << num[idx].LaTeX();
        else ss << num[idx].LaTeX_sum();
      }
    }
    else{
      ss << "\\dfrac{" ;
      for (unsigned int idx=0;idx<num.size();idx++) {
        if (idx==0) ss << num[idx].LaTeX();
        else ss << num[idx].LaTeX_sum();
      }
    ss<< "}{" << den.LaTeX() << "}";
    }
    return ss.str();
    }  
  };

std::vector<expr_product> parse_eqn(const std::string& eq){
  std::vector<expr_product> ans;
  std::string substring;
  for (std::string::const_iterator itr=eq.begin();itr!=eq.end();itr++){
    if (*itr=='-'||*itr=='+') {ans.push_back(substring); substring=*itr;}
    else substring.push_back(*itr);
  }
  ans.push_back(substring);
  std::vector<expr_product> sane_ans;
  for (int idx=0;idx<ans.size();idx++) if (!ans[idx].is_zero()) sane_ans.push_back(ans[idx]);
  return sane_ans;
  }

std::vector<expr_product> equalize_expr(const std::string& formula){
  std::vector<expr_product> ans;
  std::string left,right;
  std::string::const_iterator cursor;
  cursor = std::find(formula.begin(),formula.end(),'=');
  if (cursor==formula.end()) return parse_eqn(formula);
  cursor = formula.begin();
  while (*cursor!='='&&cursor!=formula.end()) left.push_back(*cursor++);
  if (*cursor=='=') cursor++;
  while (cursor!=formula.end()) right.push_back(*cursor++);
  std::vector<expr_product> right_part;
  ans = parse_eqn(left);
  right_part = parse_eqn(right);
  for (int idx=0;idx<right_part.size();idx++) ans.push_back(-right_part[idx]);
  return ans;
}

std::string equalize(const std::string& formula){
  std::vector<expr_product> parts = equalize_expr(formula);
  std::string ans;
  for (int idx=0;idx<parts.size();idx++) ans+=" "+parts[idx].str();
  return ans;
  }

expr_fraction solve(const std::string& formula,const std::string& variable){
  //expr_fraction ans;
  std::vector<expr_product> equation = equalize_expr(formula);
  //std::cout << std::endl;
  //std::cout << "Flattened equation:" ;
  //for (int idx=0;idx<equation.size();idx++) std::cout << equation[idx];
  //std::cout << std::endl;
  
  int index=equation.size();
  for (int idx=0;idx<equation.size();idx++)
    if (equation[idx].GetA()==variable||equation[idx].GetB()==variable) {
      index=idx; 
      break;}
  //std::cout << "Term with variable(" << variable << ") is " << index << std::endl;
  //std::cout << equation[index] << std::endl;
  if (equation[index].GetA()!=variable)  equation[index].exchange();
  //std::cout << "After Permutation:" << equation[index] << std::endl;
  std::vector<expr_product> num;
  expr_product den(equation[index].GetB());
  //std::cout << "Denominator = " << den << std::endl;
  //std::cout << "Numerator =";
  for (int idx=0;idx<equation.size();idx++){
    if (idx!=index) {
      //if (equation[index].is_negative())  std::cout << equation[idx] << " " ;
      //else                                std::cout << -equation[idx] << " " ;
      if (equation[index].is_negative())  num.push_back(equation[idx]);
      else                                num.push_back(-equation[idx]);
    } // end if
  }// end for
  return expr_fraction(num,den);
}

void echo_LU(const std::string& mode,const int& n){
  matrix<std::string> A,L,u,p;
  bool reverse=false;
  if (mode=="Ul"||mode=="uL") reverse=true;
  A = build_lu_matrix(n,"a");
  L = build_lu_matrix(n,std::string(1,mode[0]));
  u = build_lu_matrix(n,std::string(1,mode[1]));
  p=L*u;
  
 for (unsigned int idy=0;idy<n;idy++){
   std::cout << "(" << std::flush;
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << A(idx,idy) << std::flush;
   std::cout << " ) = ( " << std::flush;
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << L(idx,idy) << std::flush;
   std::cout << " ) * ( " << std::flush;
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << u(idx,idy) << std::flush;
   std::cout << " ) " << std::endl;
 }
  if (reverse){ A=A.reverse(); p=p.reverse();} 
  for (unsigned int idx=0;idx<n;idx++) 
  for (unsigned int idy=0;idy<n;idy++)
  std::cout << A(idx,idy) << " = "  << p(idx,idy) << std::endl;
  return;
  }

void LaTeX_LU(const std::string& mode,const int& n){
  matrix<std::string> A,L,U,p,LU;
  bool reverse=false;
  std::string lmode,umode;
  lmode+=mode[0]; lmode+=" t"; umode+=mode[1]; umode+=" t";
  if (mode=="Ul"||mode=="uL") reverse=true;
  A = build_lu_matrix(n,"a t");
  L = build_lu_matrix(n,lmode);
  U = build_lu_matrix(n,umode);
  p = L*U;
  LU = build_lu_matrix(n,mode+" t");
  
  std::cout << "\\subsection{Caso de matrices de $"<< n<<"\\times "<< n <<"$ }" << std::endl;
  
  std::cout << "$$ " << A.LaTeX() << " = " << L.LaTeX() << " \\cdot " << U.LaTeX() << " $$" << std::endl;
  std::cout << "$$ " << LU.LaTeX() << " $$" << std::endl;
  
  std::cout << "\\begin{multicols}{2}" << std::endl;
  if (reverse){ A=A.reverse(); p=p.reverse();LU=LU.reverse();}
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++)
    std::cout << "$$ " << A(idx,idy) << " = " << translate(p(idx,idy),'*',' ') << " $$" << std::endl;
  std::cout << "\\vfill\\null" << std::endl << "\\columnbreak" << std::endl;
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++) {
    std::cout << "$$ " << LU(idx,idy) << " = ";
    std::cout << solve(A(idx,idy)+" = "+p(idx,idy),LU(idx,idy)).LaTeX();
    std::cout << " $$" << std::endl;
  } std::cout << "\\end{multicols}" <<std::endl;
  return;
}

void LaTeX_LU_solve(const std::string& mode,const int& n){
  matrix<std::string> L,x,y,p;
  bool reverse=false;
  if (mode=="U"||mode=="u") reverse=true;
  L = build_lu_matrix(n,mode+" t");
  x = build_row_matrix(n,"x t");
  y = build_row_matrix(n,"y t");
  std::cout << "\\subsection{Caso de matrices de $"<< n<<"\\times "<< n <<"$ }" << std::endl;
  std::cout << "$$ " << y.LaTeX() << " = " << L.LaTeX() << " \\cdot " << x.LaTeX() << " $$" << std::endl;
  p = L*x;
  std::cout << "\\begin{multicols}{2}" << std::endl;
  if (reverse){ y=y.reverse(); p=p.reverse();x=x.reverse();}
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ " << y(0,idx) << " = " << translate(p(0,idx),'*',' ') << " $$" << std::endl;
  std::cout << "\\vfill\\null" << std::endl << "\\columnbreak" << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ "  << x(0,idx) << " = "<< solve(y(0,idx)+" = "+p(0,idx),x(0,idx)).LaTeX() << " $$" << std::endl;
  std::cout << "\\end{multicols}" <<std::endl;
  return;
}

void echo_Lu(const int& n){ echo_LU("Lu",n); return;}
void echo_lU(const int& n){ echo_LU("lU",n); return;}
void echo_uL(const int& n){ echo_LU("uL",n); return;}
void echo_Ul(const int& n){ echo_LU("Ul",n); return;}

void LaTeX_Lu(const int& n){ LaTeX_LU("Lu",n); return;}
void LaTeX_lU(const int& n){ LaTeX_LU("lU",n); return;}
void LaTeX_uL(const int& n){ LaTeX_LU("uL",n); return;}
void LaTeX_Ul(const int& n){ LaTeX_LU("Ul",n); return;}

void LaTeX_solve_L(const int& n){ LaTeX_LU_solve("L",n); return;}
void LaTeX_solve_l(const int& n){ LaTeX_LU_solve("l",n); return;}
void LaTeX_solve_U(const int& n){ LaTeX_LU_solve("U",n); return;}
void LaTeX_solve_u(const int& n){ LaTeX_LU_solve("u",n); return;}

void tests(){
  hbar();
  echo_Lu(4);
  hbar();
  echo_lU(4);
  hbar();
  echo_Ul(4);
  hbar();
  echo_uL(4);
  hbar();
  expr_product test;
  test.parse_str(std::string(" - xx * yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string("-xx*yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string(" - xx*yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string("- xx * yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  hbar();
  test.parse_str(std::string("  xx * yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string(" + xx * yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string("+xx*yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string(" + xx*yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string("+ xx * yy"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  hbar();
  test.parse_str(std::string("  xx "));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string("  -xx   "));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  test.parse_str(std::string(" xx * 1"));
  std::cout << "Testing expr_product:" <<test << std::endl; 
  hbar();
  std::cout << "Equalize x*y + a*b = z : " << equalize("x*y + a*b = z")  << std::endl;
  std::cout << "Equalize t*z = x*y + a*b : " << equalize("t*z = x*y + a*b")  << std::endl;
  hbar();
  expr_fraction solve_test;
  std::cout << "Solve x in a*x + b*y + c*z = 0 : " << solve("a*x + b*y + c*z = 0","x").str() << std::endl;
  std::cout << "Solve x in -a*x + b*y + c*z = 0 : " << solve("-a*x + b*y + c*z = 0","x").str() << std::endl;
  std::cout << "Solve x in a*x - b*y - c*z = 0 : " << solve("a*x - b*y - c*z = 0","x").str() << std::endl;
  std::cout << "Solve x in t*z = x*y + a*b : " << solve("t*z = x*y + a*b","x").str()  << std::endl;
  std::cout << "Solve x in t*x = z*y + a*b : " << solve("t*x = z*y + a*b","x").str()  << std::endl;
  std::cout << "Solve x in t*z = y*x + a*b : " << solve("t*z = y*x + a*b","x").str()  << std::endl;
  std::cout << "Solve x in x*t = z*y + a*b : " << solve("x*t = z*y + a*b","x").str()  << std::endl;
  hbar();
  return;
  }

void matrix_tests(){
  std::cout << build_lu_matrix(4,"L LaTeX").LaTeX() << "\\cdot " ;
  std::cout << build_lu_matrix(4,"u LaTeX").LaTeX() << " = " ;
  std::cout << translate((build_lu_matrix(4,"L LaTeX")*build_lu_matrix(4,"u LaTeX")).LaTeX(),'*',' ');
  std::cout << std::endl;
  return;
}

int main(){
  std::vector<std::string> LU_modes,L_modes; 
  LU_modes.push_back("Lu");  LU_modes.push_back("lU");  LU_modes.push_back("uL");  LU_modes.push_back("Ul");
  L_modes.push_back("L"); L_modes.push_back("l"); L_modes.push_back("U"); L_modes.push_back("u");
  std::cout << "\\documentclass[10pt,a4paper,dvipdfmx]{article}" << std::endl;
  std::cout << "\\usepackage[utf8]{inputenc}" << std::endl;
  std::cout << "\\usepackage[english]{babel}" << std::endl;
  std::cout << "\\usepackage{amsmath}" << std::endl;
  std::cout << "\\usepackage{amsfonts}" << std::endl;
  std::cout << "\\usepackage{amssymb}" << std::endl;
  std::cout << "\\usepackage{multicol}" << std::endl;
  std::cout << "\\usepackage[hidelinks,dvipdfm]{hyperref}" << std::endl;
  std::cout << "\\usepackage{graphicx}" << std::endl;
  std::cout << "\\usepackage{kpfonts}" << std::endl;
  std::cout << "\\usepackage[left=1.6cm,right=1.4cm,top=1.5cm,bottom=1.5cm]{geometry}" << std::endl;
  std::cout << "\\title{M\\\'etodos para resoluci\\\'on de sistemas de ecuaciones lineales $LU$}" << std::endl;
  std::cout << "\\author{Ricardo-Francisco Luis Mart\\\'inez}" << std::endl;
  std::cout << "\\begin{document}" << std::endl ;
  std::cout << "\\maketitle" << std::endl;
  std::cout << "\\tableofcontents" << std::endl;
  std::cout << "\\newpage" << std::endl;
  for (int idy=0;idy<4;idy++){
    std::cout << "\\section{Metodo " << LU_modes[idy] << "}" << std::endl;
    for (int idx=2;idx<7;idx++) LaTeX_LU(LU_modes[idy],idx);
  }
  for (int idy=0;idy<4;idy++){
    std::cout << "\\section{Resoluci\\\'on de matrices triangulares "<< L_modes[idy]<<"}"<< std::endl;
    for (int idx=2;idx<7;idx++) LaTeX_LU_solve(L_modes[idy],idx);
  }
  std::cout << "\\end{document}" << std::endl;
  
  return 0;
  }


/*
void LaTeX_Lu(const int& n){
  matrix<std::string> A,L,U,p,LU;
  A = build_lu_matrix(n,"a t");
  L = build_lu_matrix(n,"L t");
  U = build_lu_matrix(n,"u t");
  p = L*U;
  LU = build_lu_matrix(n,"Lu t");

  std::cout << "$$ " << A.LaTeX() << " = " << L.LaTeX() << " \\cdot " << U.LaTeX() << " $$" << std::endl;
  std::cout << "$$ " << LU.LaTeX() << " $$" << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++){
    std::cout << "$$ " << A(idx,idy) << " = " << translate(p(idx,idy),'*',' ') << " $$" << std::endl;
  }
  std::cout << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++) {
    std::cout << "$$ " << LU(idx,idy) << " = ";
    std::cout << solve(A(idx,idy)+" = "+p(idx,idy),LU(idx,idy)).LaTeX();
    std::cout << " $$" << std::endl;
  } std::cout << std::endl;
  return;
}

void LaTeX_lU(const int& n){
  matrix<std::string> A,L,U,p,LU;
  A = build_lu_matrix(n,"a t");
  L = build_lu_matrix(n,"l t");
  U = build_lu_matrix(n,"U t");
  p = L*U;
  LU = build_lu_matrix(n,"lU t");

  std::cout << "$$ " << A.LaTeX() << " = " << L.LaTeX() << " \\cdot " << U.LaTeX() << " $$" << std::endl;
  std::cout << "$$ " << LU.LaTeX() << " $$" << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++){
    std::cout << "$$ " << A(idx,idy) << " = " << translate(p(idx,idy),'*',' ') << " $$" << std::endl;
  }
  std::cout << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++) {
    std::cout << "$$ " << LU(idx,idy) << " = ";
    std::cout << solve(A(idx,idy)+" = "+p(idx,idy),LU(idx,idy)).LaTeX();
    std::cout << " $$" << std::endl;
  } std::cout << std::endl;
  return;
}

void LaTeX_uL(const int& n){
  matrix<std::string> A,L,U,p,LU;
  A = build_lu_matrix(n,"a t");
  L = build_lu_matrix(n,"L t");
  U = build_lu_matrix(n,"u t");
  p = U*L;
  LU = build_lu_matrix(n,"Lu t");

  std::cout << "$$ " << A.LaTeX() << " = " << U.LaTeX() << " \\cdot " << L.LaTeX() << " $$" << std::endl;
  std::cout << "$$ " << LU.LaTeX() << " $$" << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++){
    std::cout << "$$ " << A(n-1-idx,n-1-idy) << " = " << translate(p(n-1-idx,n-1-idy),'*',' ') << " $$" << std::endl;
  }
  std::cout << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++) {
    std::cout << "$$ " << LU(n-1-idx,n-1-idy) << " = ";
    std::cout << solve(A(n-1-idx,n-1-idy)+" = "+p(n-1-idx,n-1-idy),LU(n-1-idx,n-1-idy)).LaTeX();
    std::cout << " $$" << std::endl;
  } std::cout << std::endl;
  return;
}

void LaTeX_Ul(const int& n){
  matrix<std::string> A,L,U,p,LU;
  A = build_lu_matrix(n,"a t");
  L = build_lu_matrix(n,"l t");
  U = build_lu_matrix(n,"U t");
  p = U*L;
  LU = build_lu_matrix(n,"lU t");

  std::cout << "$$ " << A.LaTeX() << " = " << U.LaTeX() << " \\cdot " << L.LaTeX() << " $$" << std::endl;
  std::cout << "$$ " << LU.LaTeX() << " $$" << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++){
    std::cout << "$$ " << A(n-1-idx,n-1-idy) << " = " << translate(p(n-1-idx,n-1-idy),'*',' ') << " $$" << std::endl;
  }
  std::cout << std::endl;
  
  for (unsigned int idx=0;idx<n;idx++) for (unsigned int idy=0;idy<n;idy++) {
    std::cout << "$$ " << LU(n-1-idx,n-1-idy) << " = ";
    std::cout << solve(A(n-1-idx,n-1-idy)+" = "+p(n-1-idx,n-1-idy),LU(n-1-idx,n-1-idy)).LaTeX();
    std::cout << " $$" << std::endl;
  } std::cout << std::endl;
  return;
}

*/

/*
void LaTeX_solve_L(const int& n){
  matrix<std::string> L,x,y,p;
  L = build_lu_matrix(n,"L t");
  x = build_row_matrix(n,"x t");
  y = build_row_matrix(n,"y t");
  std::cout << "$$ " << y.LaTeX() << " = " << L.LaTeX() << " \\cdot " << x.LaTeX() << " $$" << std::endl;
  p = L*x;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ " << y(0,idx) << " = " << translate(p(0,idx),'*',' ') << " $$" << std::endl;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ "  << x(0,idx) << " = "<< solve(y(0,idx)+" = "+p(0,idx),x(0,idx)).LaTeX() << " $$" << std::endl;
  std::cout << std::endl;
  return;
  }
  
void LaTeX_solve_l(const int& n){
  matrix<std::string> L,x,y,p;
  L = build_lu_matrix(n,"l t");
  x = build_row_matrix(n,"x t");
  y = build_row_matrix(n,"y t");
  std::cout << "$$ " << y.LaTeX() << " = " << L.LaTeX() << " \\cdot " << x.LaTeX() << " $$" << std::endl;
  p = L*x;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ " << y(0,idx) << " = " << translate(p(0,idx),'*',' ') << " $$" << std::endl;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ "  << x(0,idx) << " = "<< solve(y(0,idx)+" = "+p(0,idx),x(0,idx)).LaTeX() << " $$" << std::endl;
  std::cout << std::endl;
  return;
  }
  
void LaTeX_solve_U(const int& n){
  matrix<std::string> L,x,y,p;
  L = build_lu_matrix(n,"U t");
  x = build_row_matrix(n,"x t");
  y = build_row_matrix(n,"y t");
  std::cout << "$$ " << y.LaTeX() << " = " << L.LaTeX() << " \\cdot " << x.LaTeX() << " $$" << std::endl;
  p = L*x;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ " << y(0,n-1-idx) << " = " << translate(p(0,n-1-idx),'*',' ') << " $$" << std::endl;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ "  << x(0,n-1-idx) << " = "<< solve(y(0,n-1-idx)+" = "+p(0,n-1-idx),x(0,n-1-idx)).LaTeX() << " $$" << std::endl;
  std::cout << std::endl;
  return;
  }
  
void LaTeX_solve_u(const int& n){
  matrix<std::string> L,x,y,p;
  L = build_lu_matrix(n,"u t");
  x = build_row_matrix(n,"x t");
  y = build_row_matrix(n,"y t");
  std::cout << "$$ " << y.LaTeX() << " = " << L.LaTeX() << " \\cdot " << x.LaTeX() << " $$" << std::endl;
  p = L*x;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ " << y(0,n-1-idx) << " = " << translate(p(0,n-1-idx),'*',' ') << " $$" << std::endl;
  std::cout << std::endl;
  for (unsigned int idx=0;idx<n;idx++)
    std::cout << "$$ "  << x(0,n-1-idx) << " = "<< solve(y(0,n-1-idx)+" = "+p(0,n-1-idx),x(0,n-1-idx)).LaTeX() << " $$" << std::endl;
  std::cout << std::endl;
  return;
  }  

*/

/*
void echo_Lu(const int& n){
  matrix<std::string> A,L,u,p;
  A = build_lu_matrix(n,"a");
  L = build_lu_matrix(n,"L");
  u = build_lu_matrix(n,"u");
  p=L*u;
 
 
  
 for (unsigned int idy=0;idy<n;idy++){
   std::cout << "(" << std::flush;
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << A(idx,idy) << std::flush;
   std::cout << " ) = ( " << std::flush;
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << L(idx,idy) << std::flush;
   std::cout << " ) * ( " << std::flush;
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << u(idx,idy) << std::flush;
   std::cout << " ) " << std::endl;
 }
 
 
  for (unsigned int idx=0;idx<n;idx++) 
  for (unsigned int idy=0;idy<n;idy++)
  std::cout << A(idx,idy) << " = "  << p(idx,idy) << std::endl;
  return;
}

void echo_lU(const int& n){
  matrix<std::string> A,L,u,p;
  L = build_lu_matrix(n,"l");
  u = build_lu_matrix(n,"U");
  p=L*u;
  
 for (unsigned int idy=0;idy<n;idy++){
   std::cout << "(";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << A(idx,idy);
   std::cout << " ) = ( ";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << L(idx,idy);
   std::cout << " ) * ( ";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << u(idx,idy);
   std::cout << " ) " << std::endl;
 }
 
 
  for (unsigned int idx=0;idx<n;idx++) 
  for (unsigned int idy=0;idy<n;idy++)
  std::cout << A(idx,idy) << " = "  << p(idx,idy) << std::endl;
  return;
}

void echo_uL(const int& n){
  matrix<std::string> A,L,u,p;
  L = build_lu_matrix(n,"L");
  u = build_lu_matrix(n,"u");
  p=u*L;
  
 for (unsigned int idy=0;idy<n;idy++){
   std::cout << "(";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << A(idx,idy);
   std::cout << " ) = ( ";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << L(idx,idy);
   std::cout << " ) * ( ";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << u(idx,idy);
   std::cout << " ) " << std::endl;
 }
 
 
  for (unsigned int idx=0;idx<n;idx++) 
  for (unsigned int idy=0;idy<n;idy++)
  std::cout << A(n-1-idx,n-1-idy) << " = "  << p(n-1-idx,n-1-idy) << std::endl;
  return;
}

void echo_Ul(const int& n){
  matrix<std::string> A,L,u,p;
  L = build_lu_matrix(n,"l");
  u = build_lu_matrix(n,"U");
  p=u*L;
  
 for (unsigned int idy=0;idy<n;idy++){
   std::cout << "(";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << A(idx,idy);
   std::cout << " ) = ( ";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << L(idx,idy);
   std::cout << " ) * ( ";
   for (unsigned int idx=0;idx<n;idx++) 
   std::cout << std::setw(4) << u(idx,idy);
   std::cout << " ) " << std::endl;
  }
 
  for (unsigned int idx=0;idx<n;idx++) 
  for (unsigned int idy=0;idy<n;idy++)
  std::cout << A(n-1-idx,n-1-idy) << " = "  << p(n-1-idx,n-1-idy) << std::endl;
  return;
}

*/


/*
 */
