#ifndef _TID_CHECKS_HPP
#define _TID_CHECKS_HPP

#include <stdint.h>

#define TID_T uint32_t

/* attemp to centralize any references to harcoded tax ids here */

#define ART_SEQ_TID 32630
#define HUMAN_TID 9606

#define isPhiX(tid) ( (tid==374840)||(tid==10847)||(tid==374840)||(tid==32630) ? 1 : 0)

inline bool isHuman(TID_T taxid)  {
   bool res = false;
   switch (taxid) {
   case 9606:
   case 63221: //neanderthal:
   case 741158: //denisovan
      res = true;
      break;
   default:
      res=false;
      break;
   }
   return res;
}

inline bool isEukId(TID_T taxid) {
   bool res = false;
   switch (taxid) {
   case 2759:
      res = true;
      break;
   default:
      res = false;
      break; 
   }
   return res;
}

inline bool isVirId(TID_T taxid) {
   bool res = false;
   switch (taxid) {
   case 10239:
      res = true;
      break;
   default:
      res = false;
      break;
   }
   return res;
}

inline bool isProkId(TID_T taxid) {
   bool res = false;
   switch (taxid) {
   case 2157: //Archaea
   case 2:
      res = true;
      break;
   default:
      res = false;
      break;
   }
   return res;
}


#endif // _TID_CHECKS_HPP
