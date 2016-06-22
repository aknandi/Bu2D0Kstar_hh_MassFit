// Include files 

#include "Base.h"
#include "CommonTools.h"

//-----------------------------------------------------------------------------
// Implementation file for class : Base
//-----------------------------------------------------------------------------

Base::Base()
  : d2kpi("d2kpi")
  , d2kk("d2kk")
  , d2pipi("d2pipi")
  , d2pik("d2pik")
  , minus("minus")
  , plus("plus")
  , both("both")
  , LL("LL")
  , DD("DD")
  , mix("mix")
  , run1("run1")
  , run2("run2")
  , all("all")
  , merge("merge")
  , separate("separate")
  , slash("/")
  , underscore("_")
  , null("NULL")
  , pi(3.141592653589793)
  , twopi(6.283185307179586)
  , pibytwo(1.570796326794897)
{

  allmodeList.push_back(d2kpi);
  allmodeList.push_back(d2kk);
  allmodeList.push_back(d2pipi);
  allmodeList.push_back(d2pik);
        
}
