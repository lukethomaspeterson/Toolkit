#include "arit-r.h"

class complex
{
   renum re;
   renum im;

public:

   renum real() const;
   renum imag() const;

   complex();
   complex(const complex& y);
   complex(renum r, renum i=0);

   complex& operator =  (const complex& y);
   complex& operator =  (const renum& y);

   complex& operator += (const complex& y);
   complex& operator += (renum y);
   complex& operator -= (const complex& y);
   complex& operator -= (renum y);
   complex& operator *= (const complex& y);
   complex& operator *= (renum y);

   complex& operator /= (const complex& y); 
   complex& operator /= (renum y); 

   void error(const char* msg) const;
};

// other functions defined as inlines

complex  operator - (const complex& x);
complex  conj(const complex& x);
complex  operator + (const complex& x, const complex& y);
complex  operator + (const complex& x, renum y);
complex  operator + (renum x, const complex& y);
complex  operator - (const complex& x, const complex& y);
complex  operator - (const complex& x, renum y);
complex  operator - (renum x, const complex& y);
complex  operator * (const complex& x, const complex& y);
complex  operator * (const complex& x, renum y);
complex  operator * (const complex& x, int y);
complex  operator * (renum x, const complex& y);

renum real(const complex& x);
renum imag(const complex& x);

// inline members

inline renum complex::real() const { return re; }
inline renum complex::imag() const { return im; }

inline complex::complex() {}
inline complex::complex(const complex& y) :re(y.real()), im(y.imag()) {}
inline complex::complex(renum r, renum i) :re(r), im(i) {}

inline complex& complex::operator = (const complex& y) 
{ 
  re = y.real(); im = y.imag(); return *this; 
} 

inline complex& complex::operator = (const renum& y) 
{ 
  re = y; im = 0; return *this; 
} 

inline complex& complex::operator += (const complex& y)
{ 
  re += y.real();  im += y.imag(); return *this; 
}

inline complex& complex::operator += (renum y)
{ 
  re += y; return *this; 
}

inline complex& complex::operator -= (const complex& y)
{ 
  re -= y.real();  im -= y.imag(); return *this; 
}

inline complex& complex::operator -= (renum y)
{ 
  re -= y; return *this; 
}

inline complex& complex::operator *= (const complex& y)
{  
  renum r = re * y.real() - im * y.imag();
  im = re * y.imag() + im * y.real(); 
  re = r; 
  return *this; 
}

inline complex& complex::operator *= (renum y)
{  
  re = re * y; im = im * y; return *this; 
}

//  functions

inline complex operator - (const complex& x)
{
  return complex(-x.real(), -x.imag());
}

inline complex conj(const complex& x)
{
  return complex(x.real(), -x.imag());
}

inline complex operator + (const complex& x, const complex& y)
{
  return complex(x.real() + y.real(), x.imag() + y.imag());
}

inline complex operator + (const complex& x, renum y)
{
  return complex(x.real() + y, x.imag());
}

inline complex operator + (renum x, const complex& y)
{
  return complex(x + y.real(), y.imag());
}

inline complex operator - (const complex& x, const complex& y)
{
  return complex(x.real() - y.real(), x.imag() - y.imag());
}

inline complex operator - (const complex& x, renum y)
{
  return complex(x.real() - y, x.imag());
}

inline complex operator - (renum x, const complex& y)
{
  return complex(x - y.real(), -y.imag());
}

inline complex operator * (const complex& x, const complex& y)
{
  return complex(x.real() * y.real() - x.imag() * y.imag(), 
                 x.real() * y.imag() + x.imag() * y.real());
}

inline complex operator * (const complex& x, renum y)
{
  return complex(x.real() * y, x.imag() * y);
}

inline complex operator * (const complex& x, int y)
{
  return complex(x.real() * y, x.imag() * y);
}

inline complex operator * (renum x, const complex& y)
{
  return complex(x * y.real(), x * y.imag());
}

inline renum real(const complex& x)
{
  return x.real();
}

inline renum imag(const complex& x)
{
  return x.imag();
}


inline complex operator / (const complex& x, const complex& y)
{
  renum d=y.real()*y.real()+y.imag()*y.imag();
  return(complex((x.real()*y.real()+x.imag()*y.imag())/d,
  (x.imag()*y.real()-x.real()*y.imag())/d));
}
inline complex operator / (const complex& x, renum r)
{
   return(complex(x.real()/r,x.imag()/r));
}
