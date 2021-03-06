//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun 14 15:33:00 2011 by ROOT version 5.31/01)
//      from the StreamerInfo in file http://root.cern.ch/files/atlasFlushed.root
//////////////////////////////////////////////////////////


#ifndef Trk__Surface_p1_h
#define Trk__Surface_p1_h
namespace Trk {
class Surface_p1;
} // end of namespace.

#include "Riostream.h"
#include <vector>

namespace Trk {
class Surface_p1 {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int m_associatedDetElementId;    //
   vector<float> m_transform;                 //

   Surface_p1();
   Surface_p1(const Surface_p1 & );
   virtual ~Surface_p1();

};
} // namespace
#endif
