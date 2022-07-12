// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIjunofsdIusersdImiaoyudIsupernovadIwenljdIsimulationdIrootdictdIlibSNsim_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNchannelNuE.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnupLS.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNsource.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnumNakazatoSrc.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNdetect.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNchannelNCC.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNpreGuoIntegFcn.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnumBurrowsSrc.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNeffectLS.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnumGarchingSrc.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNGarchingIntegFcn.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNcncLS.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnccLS.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnumJapanSrc.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnumPreGuoSrc.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNanaSrc.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNchannelBCC.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNchannelCNC.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnueLS.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNpreKatoIntegFcn.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNNakazatoIntegFcn.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNchannelIBD.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNchannelNuP.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNchannels.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNbccLS.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNnumPreKatoSrc.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNJapanIntegFcn.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNibdLS.hh"
#include "/junofs/users/miaoyu/supernova/wenlj/simulation/include/SNBurrowsIntegFcn.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_SNeffectLS(void *p);
   static void deleteArray_SNeffectLS(void *p);
   static void destruct_SNeffectLS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNeffectLS*)
   {
      ::SNeffectLS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNeffectLS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNeffectLS", ::SNeffectLS::Class_Version(), "SNeffectLS.hh", 13,
                  typeid(::SNeffectLS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNeffectLS::Dictionary, isa_proxy, 4,
                  sizeof(::SNeffectLS) );
      instance.SetDelete(&delete_SNeffectLS);
      instance.SetDeleteArray(&deleteArray_SNeffectLS);
      instance.SetDestructor(&destruct_SNeffectLS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNeffectLS*)
   {
      return GenerateInitInstanceLocal((::SNeffectLS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNeffectLS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_SNchannels(void *p);
   static void deleteArray_SNchannels(void *p);
   static void destruct_SNchannels(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNchannels*)
   {
      ::SNchannels *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNchannels >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNchannels", ::SNchannels::Class_Version(), "SNchannels.hh", 6,
                  typeid(::SNchannels), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNchannels::Dictionary, isa_proxy, 4,
                  sizeof(::SNchannels) );
      instance.SetDelete(&delete_SNchannels);
      instance.SetDeleteArray(&deleteArray_SNchannels);
      instance.SetDestructor(&destruct_SNchannels);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNchannels*)
   {
      return GenerateInitInstanceLocal((::SNchannels*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNchannels*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnueLS(void *p = 0);
   static void *newArray_SNnueLS(Long_t size, void *p);
   static void delete_SNnueLS(void *p);
   static void deleteArray_SNnueLS(void *p);
   static void destruct_SNnueLS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnueLS*)
   {
      ::SNnueLS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnueLS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnueLS", ::SNnueLS::Class_Version(), "SNnueLS.hh", 14,
                  typeid(::SNnueLS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnueLS::Dictionary, isa_proxy, 4,
                  sizeof(::SNnueLS) );
      instance.SetNew(&new_SNnueLS);
      instance.SetNewArray(&newArray_SNnueLS);
      instance.SetDelete(&delete_SNnueLS);
      instance.SetDeleteArray(&deleteArray_SNnueLS);
      instance.SetDestructor(&destruct_SNnueLS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnueLS*)
   {
      return GenerateInitInstanceLocal((::SNnueLS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnueLS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNchannelNuE(void *p = 0);
   static void *newArray_SNchannelNuE(Long_t size, void *p);
   static void delete_SNchannelNuE(void *p);
   static void deleteArray_SNchannelNuE(void *p);
   static void destruct_SNchannelNuE(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNchannelNuE*)
   {
      ::SNchannelNuE *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNchannelNuE >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNchannelNuE", ::SNchannelNuE::Class_Version(), "SNchannelNuE.hh", 6,
                  typeid(::SNchannelNuE), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNchannelNuE::Dictionary, isa_proxy, 4,
                  sizeof(::SNchannelNuE) );
      instance.SetNew(&new_SNchannelNuE);
      instance.SetNewArray(&newArray_SNchannelNuE);
      instance.SetDelete(&delete_SNchannelNuE);
      instance.SetDeleteArray(&deleteArray_SNchannelNuE);
      instance.SetDestructor(&destruct_SNchannelNuE);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNchannelNuE*)
   {
      return GenerateInitInstanceLocal((::SNchannelNuE*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNchannelNuE*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnupLS(void *p = 0);
   static void *newArray_SNnupLS(Long_t size, void *p);
   static void delete_SNnupLS(void *p);
   static void deleteArray_SNnupLS(void *p);
   static void destruct_SNnupLS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnupLS*)
   {
      ::SNnupLS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnupLS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnupLS", ::SNnupLS::Class_Version(), "SNnupLS.hh", 14,
                  typeid(::SNnupLS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnupLS::Dictionary, isa_proxy, 4,
                  sizeof(::SNnupLS) );
      instance.SetNew(&new_SNnupLS);
      instance.SetNewArray(&newArray_SNnupLS);
      instance.SetDelete(&delete_SNnupLS);
      instance.SetDeleteArray(&deleteArray_SNnupLS);
      instance.SetDestructor(&destruct_SNnupLS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnupLS*)
   {
      return GenerateInitInstanceLocal((::SNnupLS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnupLS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_SNsource(void *p);
   static void deleteArray_SNsource(void *p);
   static void destruct_SNsource(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNsource*)
   {
      ::SNsource *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNsource >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNsource", ::SNsource::Class_Version(), "SNsource.hh", 12,
                  typeid(::SNsource), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNsource::Dictionary, isa_proxy, 4,
                  sizeof(::SNsource) );
      instance.SetDelete(&delete_SNsource);
      instance.SetDeleteArray(&deleteArray_SNsource);
      instance.SetDestructor(&destruct_SNsource);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNsource*)
   {
      return GenerateInitInstanceLocal((::SNsource*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNsource*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_SNdetect(void *p);
   static void deleteArray_SNdetect(void *p);
   static void destruct_SNdetect(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNdetect*)
   {
      ::SNdetect *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNdetect >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNdetect", ::SNdetect::Class_Version(), "SNdetect.hh", 15,
                  typeid(::SNdetect), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNdetect::Dictionary, isa_proxy, 4,
                  sizeof(::SNdetect) );
      instance.SetDelete(&delete_SNdetect);
      instance.SetDeleteArray(&deleteArray_SNdetect);
      instance.SetDestructor(&destruct_SNdetect);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNdetect*)
   {
      return GenerateInitInstanceLocal((::SNdetect*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNdetect*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnccLS(void *p = 0);
   static void *newArray_SNnccLS(Long_t size, void *p);
   static void delete_SNnccLS(void *p);
   static void deleteArray_SNnccLS(void *p);
   static void destruct_SNnccLS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnccLS*)
   {
      ::SNnccLS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnccLS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnccLS", ::SNnccLS::Class_Version(), "SNnccLS.hh", 14,
                  typeid(::SNnccLS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnccLS::Dictionary, isa_proxy, 4,
                  sizeof(::SNnccLS) );
      instance.SetNew(&new_SNnccLS);
      instance.SetNewArray(&newArray_SNnccLS);
      instance.SetDelete(&delete_SNnccLS);
      instance.SetDeleteArray(&deleteArray_SNnccLS);
      instance.SetDestructor(&destruct_SNnccLS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnccLS*)
   {
      return GenerateInitInstanceLocal((::SNnccLS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnccLS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNchannelNCC(void *p = 0);
   static void *newArray_SNchannelNCC(Long_t size, void *p);
   static void delete_SNchannelNCC(void *p);
   static void deleteArray_SNchannelNCC(void *p);
   static void destruct_SNchannelNCC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNchannelNCC*)
   {
      ::SNchannelNCC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNchannelNCC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNchannelNCC", ::SNchannelNCC::Class_Version(), "SNchannelNCC.hh", 7,
                  typeid(::SNchannelNCC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNchannelNCC::Dictionary, isa_proxy, 4,
                  sizeof(::SNchannelNCC) );
      instance.SetNew(&new_SNchannelNCC);
      instance.SetNewArray(&newArray_SNchannelNCC);
      instance.SetDelete(&delete_SNchannelNCC);
      instance.SetDeleteArray(&deleteArray_SNchannelNCC);
      instance.SetDestructor(&destruct_SNchannelNCC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNchannelNCC*)
   {
      return GenerateInitInstanceLocal((::SNchannelNCC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNchannelNCC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *SNpreGuoIntegFcn_Dictionary();
   static void SNpreGuoIntegFcn_TClassManip(TClass*);
   static void *new_SNpreGuoIntegFcn(void *p = 0);
   static void *newArray_SNpreGuoIntegFcn(Long_t size, void *p);
   static void delete_SNpreGuoIntegFcn(void *p);
   static void deleteArray_SNpreGuoIntegFcn(void *p);
   static void destruct_SNpreGuoIntegFcn(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNpreGuoIntegFcn*)
   {
      ::SNpreGuoIntegFcn *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SNpreGuoIntegFcn));
      static ::ROOT::TGenericClassInfo 
         instance("SNpreGuoIntegFcn", "SNpreGuoIntegFcn.hh", 6,
                  typeid(::SNpreGuoIntegFcn), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SNpreGuoIntegFcn_Dictionary, isa_proxy, 4,
                  sizeof(::SNpreGuoIntegFcn) );
      instance.SetNew(&new_SNpreGuoIntegFcn);
      instance.SetNewArray(&newArray_SNpreGuoIntegFcn);
      instance.SetDelete(&delete_SNpreGuoIntegFcn);
      instance.SetDeleteArray(&deleteArray_SNpreGuoIntegFcn);
      instance.SetDestructor(&destruct_SNpreGuoIntegFcn);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNpreGuoIntegFcn*)
   {
      return GenerateInitInstanceLocal((::SNpreGuoIntegFcn*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNpreGuoIntegFcn*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SNpreGuoIntegFcn_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SNpreGuoIntegFcn*)0x0)->GetClass();
      SNpreGuoIntegFcn_TClassManip(theClass);
   return theClass;
   }

   static void SNpreGuoIntegFcn_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnumBurrowsSrc(void *p = 0);
   static void *newArray_SNnumBurrowsSrc(Long_t size, void *p);
   static void delete_SNnumBurrowsSrc(void *p);
   static void deleteArray_SNnumBurrowsSrc(void *p);
   static void destruct_SNnumBurrowsSrc(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnumBurrowsSrc*)
   {
      ::SNnumBurrowsSrc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnumBurrowsSrc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnumBurrowsSrc", ::SNnumBurrowsSrc::Class_Version(), "SNnumBurrowsSrc.hh", 11,
                  typeid(::SNnumBurrowsSrc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnumBurrowsSrc::Dictionary, isa_proxy, 4,
                  sizeof(::SNnumBurrowsSrc) );
      instance.SetNew(&new_SNnumBurrowsSrc);
      instance.SetNewArray(&newArray_SNnumBurrowsSrc);
      instance.SetDelete(&delete_SNnumBurrowsSrc);
      instance.SetDeleteArray(&deleteArray_SNnumBurrowsSrc);
      instance.SetDestructor(&destruct_SNnumBurrowsSrc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnumBurrowsSrc*)
   {
      return GenerateInitInstanceLocal((::SNnumBurrowsSrc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnumBurrowsSrc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnumGarchingSrc(void *p = 0);
   static void *newArray_SNnumGarchingSrc(Long_t size, void *p);
   static void delete_SNnumGarchingSrc(void *p);
   static void deleteArray_SNnumGarchingSrc(void *p);
   static void destruct_SNnumGarchingSrc(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnumGarchingSrc*)
   {
      ::SNnumGarchingSrc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnumGarchingSrc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnumGarchingSrc", ::SNnumGarchingSrc::Class_Version(), "SNnumGarchingSrc.hh", 11,
                  typeid(::SNnumGarchingSrc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnumGarchingSrc::Dictionary, isa_proxy, 4,
                  sizeof(::SNnumGarchingSrc) );
      instance.SetNew(&new_SNnumGarchingSrc);
      instance.SetNewArray(&newArray_SNnumGarchingSrc);
      instance.SetDelete(&delete_SNnumGarchingSrc);
      instance.SetDeleteArray(&deleteArray_SNnumGarchingSrc);
      instance.SetDestructor(&destruct_SNnumGarchingSrc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnumGarchingSrc*)
   {
      return GenerateInitInstanceLocal((::SNnumGarchingSrc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnumGarchingSrc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *SNGarchingIntegFcn_Dictionary();
   static void SNGarchingIntegFcn_TClassManip(TClass*);
   static void *new_SNGarchingIntegFcn(void *p = 0);
   static void *newArray_SNGarchingIntegFcn(Long_t size, void *p);
   static void delete_SNGarchingIntegFcn(void *p);
   static void deleteArray_SNGarchingIntegFcn(void *p);
   static void destruct_SNGarchingIntegFcn(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNGarchingIntegFcn*)
   {
      ::SNGarchingIntegFcn *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SNGarchingIntegFcn));
      static ::ROOT::TGenericClassInfo 
         instance("SNGarchingIntegFcn", "SNGarchingIntegFcn.hh", 7,
                  typeid(::SNGarchingIntegFcn), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SNGarchingIntegFcn_Dictionary, isa_proxy, 4,
                  sizeof(::SNGarchingIntegFcn) );
      instance.SetNew(&new_SNGarchingIntegFcn);
      instance.SetNewArray(&newArray_SNGarchingIntegFcn);
      instance.SetDelete(&delete_SNGarchingIntegFcn);
      instance.SetDeleteArray(&deleteArray_SNGarchingIntegFcn);
      instance.SetDestructor(&destruct_SNGarchingIntegFcn);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNGarchingIntegFcn*)
   {
      return GenerateInitInstanceLocal((::SNGarchingIntegFcn*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNGarchingIntegFcn*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SNGarchingIntegFcn_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SNGarchingIntegFcn*)0x0)->GetClass();
      SNGarchingIntegFcn_TClassManip(theClass);
   return theClass;
   }

   static void SNGarchingIntegFcn_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_SNcncLS(void *p = 0);
   static void *newArray_SNcncLS(Long_t size, void *p);
   static void delete_SNcncLS(void *p);
   static void deleteArray_SNcncLS(void *p);
   static void destruct_SNcncLS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNcncLS*)
   {
      ::SNcncLS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNcncLS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNcncLS", ::SNcncLS::Class_Version(), "SNcncLS.hh", 14,
                  typeid(::SNcncLS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNcncLS::Dictionary, isa_proxy, 4,
                  sizeof(::SNcncLS) );
      instance.SetNew(&new_SNcncLS);
      instance.SetNewArray(&newArray_SNcncLS);
      instance.SetDelete(&delete_SNcncLS);
      instance.SetDeleteArray(&deleteArray_SNcncLS);
      instance.SetDestructor(&destruct_SNcncLS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNcncLS*)
   {
      return GenerateInitInstanceLocal((::SNcncLS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNcncLS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnumJapanSrc(void *p = 0);
   static void *newArray_SNnumJapanSrc(Long_t size, void *p);
   static void delete_SNnumJapanSrc(void *p);
   static void deleteArray_SNnumJapanSrc(void *p);
   static void destruct_SNnumJapanSrc(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnumJapanSrc*)
   {
      ::SNnumJapanSrc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnumJapanSrc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnumJapanSrc", ::SNnumJapanSrc::Class_Version(), "SNnumJapanSrc.hh", 10,
                  typeid(::SNnumJapanSrc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnumJapanSrc::Dictionary, isa_proxy, 4,
                  sizeof(::SNnumJapanSrc) );
      instance.SetNew(&new_SNnumJapanSrc);
      instance.SetNewArray(&newArray_SNnumJapanSrc);
      instance.SetDelete(&delete_SNnumJapanSrc);
      instance.SetDeleteArray(&deleteArray_SNnumJapanSrc);
      instance.SetDestructor(&destruct_SNnumJapanSrc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnumJapanSrc*)
   {
      return GenerateInitInstanceLocal((::SNnumJapanSrc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnumJapanSrc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnumPreGuoSrc(void *p = 0);
   static void *newArray_SNnumPreGuoSrc(Long_t size, void *p);
   static void delete_SNnumPreGuoSrc(void *p);
   static void deleteArray_SNnumPreGuoSrc(void *p);
   static void destruct_SNnumPreGuoSrc(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnumPreGuoSrc*)
   {
      ::SNnumPreGuoSrc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnumPreGuoSrc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnumPreGuoSrc", ::SNnumPreGuoSrc::Class_Version(), "SNnumPreGuoSrc.hh", 10,
                  typeid(::SNnumPreGuoSrc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnumPreGuoSrc::Dictionary, isa_proxy, 4,
                  sizeof(::SNnumPreGuoSrc) );
      instance.SetNew(&new_SNnumPreGuoSrc);
      instance.SetNewArray(&newArray_SNnumPreGuoSrc);
      instance.SetDelete(&delete_SNnumPreGuoSrc);
      instance.SetDeleteArray(&deleteArray_SNnumPreGuoSrc);
      instance.SetDestructor(&destruct_SNnumPreGuoSrc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnumPreGuoSrc*)
   {
      return GenerateInitInstanceLocal((::SNnumPreGuoSrc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnumPreGuoSrc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNanaSrc(void *p = 0);
   static void *newArray_SNanaSrc(Long_t size, void *p);
   static void delete_SNanaSrc(void *p);
   static void deleteArray_SNanaSrc(void *p);
   static void destruct_SNanaSrc(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNanaSrc*)
   {
      ::SNanaSrc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNanaSrc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNanaSrc", ::SNanaSrc::Class_Version(), "SNanaSrc.hh", 12,
                  typeid(::SNanaSrc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNanaSrc::Dictionary, isa_proxy, 4,
                  sizeof(::SNanaSrc) );
      instance.SetNew(&new_SNanaSrc);
      instance.SetNewArray(&newArray_SNanaSrc);
      instance.SetDelete(&delete_SNanaSrc);
      instance.SetDeleteArray(&deleteArray_SNanaSrc);
      instance.SetDestructor(&destruct_SNanaSrc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNanaSrc*)
   {
      return GenerateInitInstanceLocal((::SNanaSrc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNanaSrc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNbccLS(void *p = 0);
   static void *newArray_SNbccLS(Long_t size, void *p);
   static void delete_SNbccLS(void *p);
   static void deleteArray_SNbccLS(void *p);
   static void destruct_SNbccLS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNbccLS*)
   {
      ::SNbccLS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNbccLS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNbccLS", ::SNbccLS::Class_Version(), "SNbccLS.hh", 14,
                  typeid(::SNbccLS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNbccLS::Dictionary, isa_proxy, 4,
                  sizeof(::SNbccLS) );
      instance.SetNew(&new_SNbccLS);
      instance.SetNewArray(&newArray_SNbccLS);
      instance.SetDelete(&delete_SNbccLS);
      instance.SetDeleteArray(&deleteArray_SNbccLS);
      instance.SetDestructor(&destruct_SNbccLS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNbccLS*)
   {
      return GenerateInitInstanceLocal((::SNbccLS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNbccLS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNchannelBCC(void *p = 0);
   static void *newArray_SNchannelBCC(Long_t size, void *p);
   static void delete_SNchannelBCC(void *p);
   static void deleteArray_SNchannelBCC(void *p);
   static void destruct_SNchannelBCC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNchannelBCC*)
   {
      ::SNchannelBCC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNchannelBCC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNchannelBCC", ::SNchannelBCC::Class_Version(), "SNchannelBCC.hh", 7,
                  typeid(::SNchannelBCC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNchannelBCC::Dictionary, isa_proxy, 4,
                  sizeof(::SNchannelBCC) );
      instance.SetNew(&new_SNchannelBCC);
      instance.SetNewArray(&newArray_SNchannelBCC);
      instance.SetDelete(&delete_SNchannelBCC);
      instance.SetDeleteArray(&deleteArray_SNchannelBCC);
      instance.SetDestructor(&destruct_SNchannelBCC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNchannelBCC*)
   {
      return GenerateInitInstanceLocal((::SNchannelBCC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNchannelBCC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNchannelCNC(void *p = 0);
   static void *newArray_SNchannelCNC(Long_t size, void *p);
   static void delete_SNchannelCNC(void *p);
   static void deleteArray_SNchannelCNC(void *p);
   static void destruct_SNchannelCNC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNchannelCNC*)
   {
      ::SNchannelCNC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNchannelCNC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNchannelCNC", ::SNchannelCNC::Class_Version(), "SNchannelCNC.hh", 7,
                  typeid(::SNchannelCNC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNchannelCNC::Dictionary, isa_proxy, 4,
                  sizeof(::SNchannelCNC) );
      instance.SetNew(&new_SNchannelCNC);
      instance.SetNewArray(&newArray_SNchannelCNC);
      instance.SetDelete(&delete_SNchannelCNC);
      instance.SetDeleteArray(&deleteArray_SNchannelCNC);
      instance.SetDestructor(&destruct_SNchannelCNC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNchannelCNC*)
   {
      return GenerateInitInstanceLocal((::SNchannelCNC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNchannelCNC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *SNpreKatoIntegFcn_Dictionary();
   static void SNpreKatoIntegFcn_TClassManip(TClass*);
   static void *new_SNpreKatoIntegFcn(void *p = 0);
   static void *newArray_SNpreKatoIntegFcn(Long_t size, void *p);
   static void delete_SNpreKatoIntegFcn(void *p);
   static void deleteArray_SNpreKatoIntegFcn(void *p);
   static void destruct_SNpreKatoIntegFcn(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNpreKatoIntegFcn*)
   {
      ::SNpreKatoIntegFcn *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SNpreKatoIntegFcn));
      static ::ROOT::TGenericClassInfo 
         instance("SNpreKatoIntegFcn", "SNpreKatoIntegFcn.hh", 8,
                  typeid(::SNpreKatoIntegFcn), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SNpreKatoIntegFcn_Dictionary, isa_proxy, 4,
                  sizeof(::SNpreKatoIntegFcn) );
      instance.SetNew(&new_SNpreKatoIntegFcn);
      instance.SetNewArray(&newArray_SNpreKatoIntegFcn);
      instance.SetDelete(&delete_SNpreKatoIntegFcn);
      instance.SetDeleteArray(&deleteArray_SNpreKatoIntegFcn);
      instance.SetDestructor(&destruct_SNpreKatoIntegFcn);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNpreKatoIntegFcn*)
   {
      return GenerateInitInstanceLocal((::SNpreKatoIntegFcn*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNpreKatoIntegFcn*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SNpreKatoIntegFcn_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SNpreKatoIntegFcn*)0x0)->GetClass();
      SNpreKatoIntegFcn_TClassManip(theClass);
   return theClass;
   }

   static void SNpreKatoIntegFcn_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_SNibdLS(void *p = 0);
   static void *newArray_SNibdLS(Long_t size, void *p);
   static void delete_SNibdLS(void *p);
   static void deleteArray_SNibdLS(void *p);
   static void destruct_SNibdLS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNibdLS*)
   {
      ::SNibdLS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNibdLS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNibdLS", ::SNibdLS::Class_Version(), "SNibdLS.hh", 14,
                  typeid(::SNibdLS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNibdLS::Dictionary, isa_proxy, 4,
                  sizeof(::SNibdLS) );
      instance.SetNew(&new_SNibdLS);
      instance.SetNewArray(&newArray_SNibdLS);
      instance.SetDelete(&delete_SNibdLS);
      instance.SetDeleteArray(&deleteArray_SNibdLS);
      instance.SetDestructor(&destruct_SNibdLS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNibdLS*)
   {
      return GenerateInitInstanceLocal((::SNibdLS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNibdLS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNchannelIBD(void *p = 0);
   static void *newArray_SNchannelIBD(Long_t size, void *p);
   static void delete_SNchannelIBD(void *p);
   static void deleteArray_SNchannelIBD(void *p);
   static void destruct_SNchannelIBD(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNchannelIBD*)
   {
      ::SNchannelIBD *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNchannelIBD >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNchannelIBD", ::SNchannelIBD::Class_Version(), "SNchannelIBD.hh", 7,
                  typeid(::SNchannelIBD), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNchannelIBD::Dictionary, isa_proxy, 4,
                  sizeof(::SNchannelIBD) );
      instance.SetNew(&new_SNchannelIBD);
      instance.SetNewArray(&newArray_SNchannelIBD);
      instance.SetDelete(&delete_SNchannelIBD);
      instance.SetDeleteArray(&deleteArray_SNchannelIBD);
      instance.SetDestructor(&destruct_SNchannelIBD);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNchannelIBD*)
   {
      return GenerateInitInstanceLocal((::SNchannelIBD*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNchannelIBD*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNchannelNuP(void *p = 0);
   static void *newArray_SNchannelNuP(Long_t size, void *p);
   static void delete_SNchannelNuP(void *p);
   static void deleteArray_SNchannelNuP(void *p);
   static void destruct_SNchannelNuP(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNchannelNuP*)
   {
      ::SNchannelNuP *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNchannelNuP >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNchannelNuP", ::SNchannelNuP::Class_Version(), "SNchannelNuP.hh", 6,
                  typeid(::SNchannelNuP), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNchannelNuP::Dictionary, isa_proxy, 4,
                  sizeof(::SNchannelNuP) );
      instance.SetNew(&new_SNchannelNuP);
      instance.SetNewArray(&newArray_SNchannelNuP);
      instance.SetDelete(&delete_SNchannelNuP);
      instance.SetDeleteArray(&deleteArray_SNchannelNuP);
      instance.SetDestructor(&destruct_SNchannelNuP);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNchannelNuP*)
   {
      return GenerateInitInstanceLocal((::SNchannelNuP*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNchannelNuP*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SNnumPreKatoSrc(void *p = 0);
   static void *newArray_SNnumPreKatoSrc(Long_t size, void *p);
   static void delete_SNnumPreKatoSrc(void *p);
   static void deleteArray_SNnumPreKatoSrc(void *p);
   static void destruct_SNnumPreKatoSrc(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNnumPreKatoSrc*)
   {
      ::SNnumPreKatoSrc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SNnumPreKatoSrc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SNnumPreKatoSrc", ::SNnumPreKatoSrc::Class_Version(), "SNnumPreKatoSrc.hh", 10,
                  typeid(::SNnumPreKatoSrc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SNnumPreKatoSrc::Dictionary, isa_proxy, 4,
                  sizeof(::SNnumPreKatoSrc) );
      instance.SetNew(&new_SNnumPreKatoSrc);
      instance.SetNewArray(&newArray_SNnumPreKatoSrc);
      instance.SetDelete(&delete_SNnumPreKatoSrc);
      instance.SetDeleteArray(&deleteArray_SNnumPreKatoSrc);
      instance.SetDestructor(&destruct_SNnumPreKatoSrc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNnumPreKatoSrc*)
   {
      return GenerateInitInstanceLocal((::SNnumPreKatoSrc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNnumPreKatoSrc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *SNJapanIntegFcn_Dictionary();
   static void SNJapanIntegFcn_TClassManip(TClass*);
   static void *new_SNJapanIntegFcn(void *p = 0);
   static void *newArray_SNJapanIntegFcn(Long_t size, void *p);
   static void delete_SNJapanIntegFcn(void *p);
   static void deleteArray_SNJapanIntegFcn(void *p);
   static void destruct_SNJapanIntegFcn(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNJapanIntegFcn*)
   {
      ::SNJapanIntegFcn *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SNJapanIntegFcn));
      static ::ROOT::TGenericClassInfo 
         instance("SNJapanIntegFcn", "SNJapanIntegFcn.hh", 6,
                  typeid(::SNJapanIntegFcn), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SNJapanIntegFcn_Dictionary, isa_proxy, 4,
                  sizeof(::SNJapanIntegFcn) );
      instance.SetNew(&new_SNJapanIntegFcn);
      instance.SetNewArray(&newArray_SNJapanIntegFcn);
      instance.SetDelete(&delete_SNJapanIntegFcn);
      instance.SetDeleteArray(&deleteArray_SNJapanIntegFcn);
      instance.SetDestructor(&destruct_SNJapanIntegFcn);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNJapanIntegFcn*)
   {
      return GenerateInitInstanceLocal((::SNJapanIntegFcn*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNJapanIntegFcn*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SNJapanIntegFcn_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SNJapanIntegFcn*)0x0)->GetClass();
      SNJapanIntegFcn_TClassManip(theClass);
   return theClass;
   }

   static void SNJapanIntegFcn_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SNBurrowsIntegFcn_Dictionary();
   static void SNBurrowsIntegFcn_TClassManip(TClass*);
   static void *new_SNBurrowsIntegFcn(void *p = 0);
   static void *newArray_SNBurrowsIntegFcn(Long_t size, void *p);
   static void delete_SNBurrowsIntegFcn(void *p);
   static void deleteArray_SNBurrowsIntegFcn(void *p);
   static void destruct_SNBurrowsIntegFcn(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SNBurrowsIntegFcn*)
   {
      ::SNBurrowsIntegFcn *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SNBurrowsIntegFcn));
      static ::ROOT::TGenericClassInfo 
         instance("SNBurrowsIntegFcn", "SNBurrowsIntegFcn.hh", 6,
                  typeid(::SNBurrowsIntegFcn), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SNBurrowsIntegFcn_Dictionary, isa_proxy, 4,
                  sizeof(::SNBurrowsIntegFcn) );
      instance.SetNew(&new_SNBurrowsIntegFcn);
      instance.SetNewArray(&newArray_SNBurrowsIntegFcn);
      instance.SetDelete(&delete_SNBurrowsIntegFcn);
      instance.SetDeleteArray(&deleteArray_SNBurrowsIntegFcn);
      instance.SetDestructor(&destruct_SNBurrowsIntegFcn);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SNBurrowsIntegFcn*)
   {
      return GenerateInitInstanceLocal((::SNBurrowsIntegFcn*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SNBurrowsIntegFcn*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SNBurrowsIntegFcn_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SNBurrowsIntegFcn*)0x0)->GetClass();
      SNBurrowsIntegFcn_TClassManip(theClass);
   return theClass;
   }

   static void SNBurrowsIntegFcn_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr SNeffectLS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNeffectLS::Class_Name()
{
   return "SNeffectLS";
}

//______________________________________________________________________________
const char *SNeffectLS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNeffectLS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNeffectLS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNeffectLS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNeffectLS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNeffectLS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNeffectLS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNeffectLS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNchannels::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNchannels::Class_Name()
{
   return "SNchannels";
}

//______________________________________________________________________________
const char *SNchannels::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannels*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNchannels::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannels*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNchannels::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannels*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNchannels::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannels*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnueLS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnueLS::Class_Name()
{
   return "SNnueLS";
}

//______________________________________________________________________________
const char *SNnueLS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnueLS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnueLS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnueLS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnueLS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnueLS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnueLS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnueLS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNchannelNuE::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNchannelNuE::Class_Name()
{
   return "SNchannelNuE";
}

//______________________________________________________________________________
const char *SNchannelNuE::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuE*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNchannelNuE::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuE*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNchannelNuE::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuE*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNchannelNuE::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuE*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnupLS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnupLS::Class_Name()
{
   return "SNnupLS";
}

//______________________________________________________________________________
const char *SNnupLS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnupLS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnupLS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnupLS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnupLS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnupLS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnupLS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnupLS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNsource::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNsource::Class_Name()
{
   return "SNsource";
}

//______________________________________________________________________________
const char *SNsource::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNsource*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNsource::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNsource*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNsource::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNsource*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNsource::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNsource*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNdetect::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNdetect::Class_Name()
{
   return "SNdetect";
}

//______________________________________________________________________________
const char *SNdetect::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNdetect*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNdetect::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNdetect*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNdetect::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNdetect*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNdetect::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNdetect*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnccLS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnccLS::Class_Name()
{
   return "SNnccLS";
}

//______________________________________________________________________________
const char *SNnccLS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnccLS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnccLS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnccLS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnccLS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnccLS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnccLS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnccLS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNchannelNCC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNchannelNCC::Class_Name()
{
   return "SNchannelNCC";
}

//______________________________________________________________________________
const char *SNchannelNCC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNCC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNchannelNCC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNCC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNchannelNCC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNCC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNchannelNCC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNCC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnumBurrowsSrc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnumBurrowsSrc::Class_Name()
{
   return "SNnumBurrowsSrc";
}

//______________________________________________________________________________
const char *SNnumBurrowsSrc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumBurrowsSrc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnumBurrowsSrc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumBurrowsSrc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnumBurrowsSrc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumBurrowsSrc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnumBurrowsSrc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumBurrowsSrc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnumGarchingSrc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnumGarchingSrc::Class_Name()
{
   return "SNnumGarchingSrc";
}

//______________________________________________________________________________
const char *SNnumGarchingSrc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumGarchingSrc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnumGarchingSrc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumGarchingSrc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnumGarchingSrc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumGarchingSrc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnumGarchingSrc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumGarchingSrc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNcncLS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNcncLS::Class_Name()
{
   return "SNcncLS";
}

//______________________________________________________________________________
const char *SNcncLS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNcncLS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNcncLS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNcncLS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNcncLS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNcncLS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNcncLS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNcncLS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnumJapanSrc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnumJapanSrc::Class_Name()
{
   return "SNnumJapanSrc";
}

//______________________________________________________________________________
const char *SNnumJapanSrc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumJapanSrc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnumJapanSrc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumJapanSrc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnumJapanSrc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumJapanSrc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnumJapanSrc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumJapanSrc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnumPreGuoSrc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnumPreGuoSrc::Class_Name()
{
   return "SNnumPreGuoSrc";
}

//______________________________________________________________________________
const char *SNnumPreGuoSrc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreGuoSrc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnumPreGuoSrc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreGuoSrc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnumPreGuoSrc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreGuoSrc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnumPreGuoSrc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreGuoSrc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNanaSrc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNanaSrc::Class_Name()
{
   return "SNanaSrc";
}

//______________________________________________________________________________
const char *SNanaSrc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNanaSrc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNanaSrc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNanaSrc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNanaSrc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNanaSrc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNanaSrc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNanaSrc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNbccLS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNbccLS::Class_Name()
{
   return "SNbccLS";
}

//______________________________________________________________________________
const char *SNbccLS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNbccLS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNbccLS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNbccLS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNbccLS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNbccLS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNbccLS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNbccLS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNchannelBCC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNchannelBCC::Class_Name()
{
   return "SNchannelBCC";
}

//______________________________________________________________________________
const char *SNchannelBCC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelBCC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNchannelBCC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelBCC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNchannelBCC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelBCC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNchannelBCC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelBCC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNchannelCNC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNchannelCNC::Class_Name()
{
   return "SNchannelCNC";
}

//______________________________________________________________________________
const char *SNchannelCNC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelCNC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNchannelCNC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelCNC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNchannelCNC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelCNC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNchannelCNC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelCNC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNibdLS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNibdLS::Class_Name()
{
   return "SNibdLS";
}

//______________________________________________________________________________
const char *SNibdLS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNibdLS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNibdLS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNibdLS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNibdLS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNibdLS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNibdLS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNibdLS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNchannelIBD::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNchannelIBD::Class_Name()
{
   return "SNchannelIBD";
}

//______________________________________________________________________________
const char *SNchannelIBD::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelIBD*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNchannelIBD::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelIBD*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNchannelIBD::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelIBD*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNchannelIBD::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelIBD*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNchannelNuP::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNchannelNuP::Class_Name()
{
   return "SNchannelNuP";
}

//______________________________________________________________________________
const char *SNchannelNuP::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuP*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNchannelNuP::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuP*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNchannelNuP::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuP*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNchannelNuP::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNchannelNuP*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SNnumPreKatoSrc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SNnumPreKatoSrc::Class_Name()
{
   return "SNnumPreKatoSrc";
}

//______________________________________________________________________________
const char *SNnumPreKatoSrc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreKatoSrc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SNnumPreKatoSrc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreKatoSrc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SNnumPreKatoSrc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreKatoSrc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SNnumPreKatoSrc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SNnumPreKatoSrc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void SNeffectLS::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNeffectLS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNeffectLS::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNeffectLS::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SNeffectLS(void *p) {
      delete ((::SNeffectLS*)p);
   }
   static void deleteArray_SNeffectLS(void *p) {
      delete [] ((::SNeffectLS*)p);
   }
   static void destruct_SNeffectLS(void *p) {
      typedef ::SNeffectLS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNeffectLS

//______________________________________________________________________________
void SNchannels::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNchannels.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNchannels::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNchannels::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SNchannels(void *p) {
      delete ((::SNchannels*)p);
   }
   static void deleteArray_SNchannels(void *p) {
      delete [] ((::SNchannels*)p);
   }
   static void destruct_SNchannels(void *p) {
      typedef ::SNchannels current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNchannels

//______________________________________________________________________________
void SNnueLS::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnueLS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnueLS::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnueLS::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnueLS(void *p) {
      return  p ? new(p) ::SNnueLS : new ::SNnueLS;
   }
   static void *newArray_SNnueLS(Long_t nElements, void *p) {
      return p ? new(p) ::SNnueLS[nElements] : new ::SNnueLS[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnueLS(void *p) {
      delete ((::SNnueLS*)p);
   }
   static void deleteArray_SNnueLS(void *p) {
      delete [] ((::SNnueLS*)p);
   }
   static void destruct_SNnueLS(void *p) {
      typedef ::SNnueLS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnueLS

//______________________________________________________________________________
void SNchannelNuE::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNchannelNuE.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNchannelNuE::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNchannelNuE::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNchannelNuE(void *p) {
      return  p ? new(p) ::SNchannelNuE : new ::SNchannelNuE;
   }
   static void *newArray_SNchannelNuE(Long_t nElements, void *p) {
      return p ? new(p) ::SNchannelNuE[nElements] : new ::SNchannelNuE[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNchannelNuE(void *p) {
      delete ((::SNchannelNuE*)p);
   }
   static void deleteArray_SNchannelNuE(void *p) {
      delete [] ((::SNchannelNuE*)p);
   }
   static void destruct_SNchannelNuE(void *p) {
      typedef ::SNchannelNuE current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNchannelNuE

//______________________________________________________________________________
void SNnupLS::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnupLS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnupLS::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnupLS::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnupLS(void *p) {
      return  p ? new(p) ::SNnupLS : new ::SNnupLS;
   }
   static void *newArray_SNnupLS(Long_t nElements, void *p) {
      return p ? new(p) ::SNnupLS[nElements] : new ::SNnupLS[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnupLS(void *p) {
      delete ((::SNnupLS*)p);
   }
   static void deleteArray_SNnupLS(void *p) {
      delete [] ((::SNnupLS*)p);
   }
   static void destruct_SNnupLS(void *p) {
      typedef ::SNnupLS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnupLS

//______________________________________________________________________________
void SNsource::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNsource.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNsource::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNsource::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SNsource(void *p) {
      delete ((::SNsource*)p);
   }
   static void deleteArray_SNsource(void *p) {
      delete [] ((::SNsource*)p);
   }
   static void destruct_SNsource(void *p) {
      typedef ::SNsource current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNsource

//______________________________________________________________________________
void SNdetect::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNdetect.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNdetect::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNdetect::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SNdetect(void *p) {
      delete ((::SNdetect*)p);
   }
   static void deleteArray_SNdetect(void *p) {
      delete [] ((::SNdetect*)p);
   }
   static void destruct_SNdetect(void *p) {
      typedef ::SNdetect current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNdetect

//______________________________________________________________________________
void SNnccLS::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnccLS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnccLS::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnccLS::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnccLS(void *p) {
      return  p ? new(p) ::SNnccLS : new ::SNnccLS;
   }
   static void *newArray_SNnccLS(Long_t nElements, void *p) {
      return p ? new(p) ::SNnccLS[nElements] : new ::SNnccLS[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnccLS(void *p) {
      delete ((::SNnccLS*)p);
   }
   static void deleteArray_SNnccLS(void *p) {
      delete [] ((::SNnccLS*)p);
   }
   static void destruct_SNnccLS(void *p) {
      typedef ::SNnccLS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnccLS

//______________________________________________________________________________
void SNchannelNCC::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNchannelNCC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNchannelNCC::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNchannelNCC::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNchannelNCC(void *p) {
      return  p ? new(p) ::SNchannelNCC : new ::SNchannelNCC;
   }
   static void *newArray_SNchannelNCC(Long_t nElements, void *p) {
      return p ? new(p) ::SNchannelNCC[nElements] : new ::SNchannelNCC[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNchannelNCC(void *p) {
      delete ((::SNchannelNCC*)p);
   }
   static void deleteArray_SNchannelNCC(void *p) {
      delete [] ((::SNchannelNCC*)p);
   }
   static void destruct_SNchannelNCC(void *p) {
      typedef ::SNchannelNCC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNchannelNCC

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNpreGuoIntegFcn(void *p) {
      return  p ? new(p) ::SNpreGuoIntegFcn : new ::SNpreGuoIntegFcn;
   }
   static void *newArray_SNpreGuoIntegFcn(Long_t nElements, void *p) {
      return p ? new(p) ::SNpreGuoIntegFcn[nElements] : new ::SNpreGuoIntegFcn[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNpreGuoIntegFcn(void *p) {
      delete ((::SNpreGuoIntegFcn*)p);
   }
   static void deleteArray_SNpreGuoIntegFcn(void *p) {
      delete [] ((::SNpreGuoIntegFcn*)p);
   }
   static void destruct_SNpreGuoIntegFcn(void *p) {
      typedef ::SNpreGuoIntegFcn current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNpreGuoIntegFcn

//______________________________________________________________________________
void SNnumBurrowsSrc::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnumBurrowsSrc.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnumBurrowsSrc::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnumBurrowsSrc::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnumBurrowsSrc(void *p) {
      return  p ? new(p) ::SNnumBurrowsSrc : new ::SNnumBurrowsSrc;
   }
   static void *newArray_SNnumBurrowsSrc(Long_t nElements, void *p) {
      return p ? new(p) ::SNnumBurrowsSrc[nElements] : new ::SNnumBurrowsSrc[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnumBurrowsSrc(void *p) {
      delete ((::SNnumBurrowsSrc*)p);
   }
   static void deleteArray_SNnumBurrowsSrc(void *p) {
      delete [] ((::SNnumBurrowsSrc*)p);
   }
   static void destruct_SNnumBurrowsSrc(void *p) {
      typedef ::SNnumBurrowsSrc current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnumBurrowsSrc

//______________________________________________________________________________
void SNnumGarchingSrc::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnumGarchingSrc.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnumGarchingSrc::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnumGarchingSrc::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnumGarchingSrc(void *p) {
      return  p ? new(p) ::SNnumGarchingSrc : new ::SNnumGarchingSrc;
   }
   static void *newArray_SNnumGarchingSrc(Long_t nElements, void *p) {
      return p ? new(p) ::SNnumGarchingSrc[nElements] : new ::SNnumGarchingSrc[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnumGarchingSrc(void *p) {
      delete ((::SNnumGarchingSrc*)p);
   }
   static void deleteArray_SNnumGarchingSrc(void *p) {
      delete [] ((::SNnumGarchingSrc*)p);
   }
   static void destruct_SNnumGarchingSrc(void *p) {
      typedef ::SNnumGarchingSrc current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnumGarchingSrc

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNGarchingIntegFcn(void *p) {
      return  p ? new(p) ::SNGarchingIntegFcn : new ::SNGarchingIntegFcn;
   }
   static void *newArray_SNGarchingIntegFcn(Long_t nElements, void *p) {
      return p ? new(p) ::SNGarchingIntegFcn[nElements] : new ::SNGarchingIntegFcn[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNGarchingIntegFcn(void *p) {
      delete ((::SNGarchingIntegFcn*)p);
   }
   static void deleteArray_SNGarchingIntegFcn(void *p) {
      delete [] ((::SNGarchingIntegFcn*)p);
   }
   static void destruct_SNGarchingIntegFcn(void *p) {
      typedef ::SNGarchingIntegFcn current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNGarchingIntegFcn

//______________________________________________________________________________
void SNcncLS::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNcncLS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNcncLS::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNcncLS::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNcncLS(void *p) {
      return  p ? new(p) ::SNcncLS : new ::SNcncLS;
   }
   static void *newArray_SNcncLS(Long_t nElements, void *p) {
      return p ? new(p) ::SNcncLS[nElements] : new ::SNcncLS[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNcncLS(void *p) {
      delete ((::SNcncLS*)p);
   }
   static void deleteArray_SNcncLS(void *p) {
      delete [] ((::SNcncLS*)p);
   }
   static void destruct_SNcncLS(void *p) {
      typedef ::SNcncLS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNcncLS

//______________________________________________________________________________
void SNnumJapanSrc::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnumJapanSrc.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnumJapanSrc::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnumJapanSrc::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnumJapanSrc(void *p) {
      return  p ? new(p) ::SNnumJapanSrc : new ::SNnumJapanSrc;
   }
   static void *newArray_SNnumJapanSrc(Long_t nElements, void *p) {
      return p ? new(p) ::SNnumJapanSrc[nElements] : new ::SNnumJapanSrc[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnumJapanSrc(void *p) {
      delete ((::SNnumJapanSrc*)p);
   }
   static void deleteArray_SNnumJapanSrc(void *p) {
      delete [] ((::SNnumJapanSrc*)p);
   }
   static void destruct_SNnumJapanSrc(void *p) {
      typedef ::SNnumJapanSrc current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnumJapanSrc

//______________________________________________________________________________
void SNnumPreGuoSrc::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnumPreGuoSrc.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnumPreGuoSrc::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnumPreGuoSrc::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnumPreGuoSrc(void *p) {
      return  p ? new(p) ::SNnumPreGuoSrc : new ::SNnumPreGuoSrc;
   }
   static void *newArray_SNnumPreGuoSrc(Long_t nElements, void *p) {
      return p ? new(p) ::SNnumPreGuoSrc[nElements] : new ::SNnumPreGuoSrc[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnumPreGuoSrc(void *p) {
      delete ((::SNnumPreGuoSrc*)p);
   }
   static void deleteArray_SNnumPreGuoSrc(void *p) {
      delete [] ((::SNnumPreGuoSrc*)p);
   }
   static void destruct_SNnumPreGuoSrc(void *p) {
      typedef ::SNnumPreGuoSrc current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnumPreGuoSrc

//______________________________________________________________________________
void SNanaSrc::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNanaSrc.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNanaSrc::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNanaSrc::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNanaSrc(void *p) {
      return  p ? new(p) ::SNanaSrc : new ::SNanaSrc;
   }
   static void *newArray_SNanaSrc(Long_t nElements, void *p) {
      return p ? new(p) ::SNanaSrc[nElements] : new ::SNanaSrc[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNanaSrc(void *p) {
      delete ((::SNanaSrc*)p);
   }
   static void deleteArray_SNanaSrc(void *p) {
      delete [] ((::SNanaSrc*)p);
   }
   static void destruct_SNanaSrc(void *p) {
      typedef ::SNanaSrc current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNanaSrc

//______________________________________________________________________________
void SNbccLS::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNbccLS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNbccLS::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNbccLS::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNbccLS(void *p) {
      return  p ? new(p) ::SNbccLS : new ::SNbccLS;
   }
   static void *newArray_SNbccLS(Long_t nElements, void *p) {
      return p ? new(p) ::SNbccLS[nElements] : new ::SNbccLS[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNbccLS(void *p) {
      delete ((::SNbccLS*)p);
   }
   static void deleteArray_SNbccLS(void *p) {
      delete [] ((::SNbccLS*)p);
   }
   static void destruct_SNbccLS(void *p) {
      typedef ::SNbccLS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNbccLS

//______________________________________________________________________________
void SNchannelBCC::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNchannelBCC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNchannelBCC::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNchannelBCC::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNchannelBCC(void *p) {
      return  p ? new(p) ::SNchannelBCC : new ::SNchannelBCC;
   }
   static void *newArray_SNchannelBCC(Long_t nElements, void *p) {
      return p ? new(p) ::SNchannelBCC[nElements] : new ::SNchannelBCC[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNchannelBCC(void *p) {
      delete ((::SNchannelBCC*)p);
   }
   static void deleteArray_SNchannelBCC(void *p) {
      delete [] ((::SNchannelBCC*)p);
   }
   static void destruct_SNchannelBCC(void *p) {
      typedef ::SNchannelBCC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNchannelBCC

//______________________________________________________________________________
void SNchannelCNC::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNchannelCNC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNchannelCNC::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNchannelCNC::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNchannelCNC(void *p) {
      return  p ? new(p) ::SNchannelCNC : new ::SNchannelCNC;
   }
   static void *newArray_SNchannelCNC(Long_t nElements, void *p) {
      return p ? new(p) ::SNchannelCNC[nElements] : new ::SNchannelCNC[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNchannelCNC(void *p) {
      delete ((::SNchannelCNC*)p);
   }
   static void deleteArray_SNchannelCNC(void *p) {
      delete [] ((::SNchannelCNC*)p);
   }
   static void destruct_SNchannelCNC(void *p) {
      typedef ::SNchannelCNC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNchannelCNC

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNpreKatoIntegFcn(void *p) {
      return  p ? new(p) ::SNpreKatoIntegFcn : new ::SNpreKatoIntegFcn;
   }
   static void *newArray_SNpreKatoIntegFcn(Long_t nElements, void *p) {
      return p ? new(p) ::SNpreKatoIntegFcn[nElements] : new ::SNpreKatoIntegFcn[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNpreKatoIntegFcn(void *p) {
      delete ((::SNpreKatoIntegFcn*)p);
   }
   static void deleteArray_SNpreKatoIntegFcn(void *p) {
      delete [] ((::SNpreKatoIntegFcn*)p);
   }
   static void destruct_SNpreKatoIntegFcn(void *p) {
      typedef ::SNpreKatoIntegFcn current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNpreKatoIntegFcn

//______________________________________________________________________________
void SNibdLS::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNibdLS.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNibdLS::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNibdLS::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNibdLS(void *p) {
      return  p ? new(p) ::SNibdLS : new ::SNibdLS;
   }
   static void *newArray_SNibdLS(Long_t nElements, void *p) {
      return p ? new(p) ::SNibdLS[nElements] : new ::SNibdLS[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNibdLS(void *p) {
      delete ((::SNibdLS*)p);
   }
   static void deleteArray_SNibdLS(void *p) {
      delete [] ((::SNibdLS*)p);
   }
   static void destruct_SNibdLS(void *p) {
      typedef ::SNibdLS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNibdLS

//______________________________________________________________________________
void SNchannelIBD::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNchannelIBD.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNchannelIBD::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNchannelIBD::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNchannelIBD(void *p) {
      return  p ? new(p) ::SNchannelIBD : new ::SNchannelIBD;
   }
   static void *newArray_SNchannelIBD(Long_t nElements, void *p) {
      return p ? new(p) ::SNchannelIBD[nElements] : new ::SNchannelIBD[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNchannelIBD(void *p) {
      delete ((::SNchannelIBD*)p);
   }
   static void deleteArray_SNchannelIBD(void *p) {
      delete [] ((::SNchannelIBD*)p);
   }
   static void destruct_SNchannelIBD(void *p) {
      typedef ::SNchannelIBD current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNchannelIBD

//______________________________________________________________________________
void SNchannelNuP::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNchannelNuP.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNchannelNuP::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNchannelNuP::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNchannelNuP(void *p) {
      return  p ? new(p) ::SNchannelNuP : new ::SNchannelNuP;
   }
   static void *newArray_SNchannelNuP(Long_t nElements, void *p) {
      return p ? new(p) ::SNchannelNuP[nElements] : new ::SNchannelNuP[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNchannelNuP(void *p) {
      delete ((::SNchannelNuP*)p);
   }
   static void deleteArray_SNchannelNuP(void *p) {
      delete [] ((::SNchannelNuP*)p);
   }
   static void destruct_SNchannelNuP(void *p) {
      typedef ::SNchannelNuP current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNchannelNuP

//______________________________________________________________________________
void SNnumPreKatoSrc::Streamer(TBuffer &R__b)
{
   // Stream an object of class SNnumPreKatoSrc.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SNnumPreKatoSrc::Class(),this);
   } else {
      R__b.WriteClassBuffer(SNnumPreKatoSrc::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNnumPreKatoSrc(void *p) {
      return  p ? new(p) ::SNnumPreKatoSrc : new ::SNnumPreKatoSrc;
   }
   static void *newArray_SNnumPreKatoSrc(Long_t nElements, void *p) {
      return p ? new(p) ::SNnumPreKatoSrc[nElements] : new ::SNnumPreKatoSrc[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNnumPreKatoSrc(void *p) {
      delete ((::SNnumPreKatoSrc*)p);
   }
   static void deleteArray_SNnumPreKatoSrc(void *p) {
      delete [] ((::SNnumPreKatoSrc*)p);
   }
   static void destruct_SNnumPreKatoSrc(void *p) {
      typedef ::SNnumPreKatoSrc current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNnumPreKatoSrc

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNJapanIntegFcn(void *p) {
      return  p ? new(p) ::SNJapanIntegFcn : new ::SNJapanIntegFcn;
   }
   static void *newArray_SNJapanIntegFcn(Long_t nElements, void *p) {
      return p ? new(p) ::SNJapanIntegFcn[nElements] : new ::SNJapanIntegFcn[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNJapanIntegFcn(void *p) {
      delete ((::SNJapanIntegFcn*)p);
   }
   static void deleteArray_SNJapanIntegFcn(void *p) {
      delete [] ((::SNJapanIntegFcn*)p);
   }
   static void destruct_SNJapanIntegFcn(void *p) {
      typedef ::SNJapanIntegFcn current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNJapanIntegFcn

namespace ROOT {
   // Wrappers around operator new
   static void *new_SNBurrowsIntegFcn(void *p) {
      return  p ? new(p) ::SNBurrowsIntegFcn : new ::SNBurrowsIntegFcn;
   }
   static void *newArray_SNBurrowsIntegFcn(Long_t nElements, void *p) {
      return p ? new(p) ::SNBurrowsIntegFcn[nElements] : new ::SNBurrowsIntegFcn[nElements];
   }
   // Wrapper around operator delete
   static void delete_SNBurrowsIntegFcn(void *p) {
      delete ((::SNBurrowsIntegFcn*)p);
   }
   static void deleteArray_SNBurrowsIntegFcn(void *p) {
      delete [] ((::SNBurrowsIntegFcn*)p);
   }
   static void destruct_SNBurrowsIntegFcn(void *p) {
      typedef ::SNBurrowsIntegFcn current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SNBurrowsIntegFcn

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 339,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETGraphmUgR_Dictionary();
   static void vectorlETGraphmUgR_TClassManip(TClass*);
   static void *new_vectorlETGraphmUgR(void *p = 0);
   static void *newArray_vectorlETGraphmUgR(Long_t size, void *p);
   static void delete_vectorlETGraphmUgR(void *p);
   static void deleteArray_vectorlETGraphmUgR(void *p);
   static void destruct_vectorlETGraphmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TGraph*>*)
   {
      vector<TGraph*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TGraph*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TGraph*>", -2, "vector", 339,
                  typeid(vector<TGraph*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETGraphmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TGraph*>) );
      instance.SetNew(&new_vectorlETGraphmUgR);
      instance.SetNewArray(&newArray_vectorlETGraphmUgR);
      instance.SetDelete(&delete_vectorlETGraphmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETGraphmUgR);
      instance.SetDestructor(&destruct_vectorlETGraphmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TGraph*> >()));

      ::ROOT::AddClassAlternate("vector<TGraph*>","std::vector<TGraph*, std::allocator<TGraph*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TGraph*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETGraphmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TGraph*>*)0x0)->GetClass();
      vectorlETGraphmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETGraphmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETGraphmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TGraph*> : new vector<TGraph*>;
   }
   static void *newArray_vectorlETGraphmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TGraph*>[nElements] : new vector<TGraph*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETGraphmUgR(void *p) {
      delete ((vector<TGraph*>*)p);
   }
   static void deleteArray_vectorlETGraphmUgR(void *p) {
      delete [] ((vector<TGraph*>*)p);
   }
   static void destruct_vectorlETGraphmUgR(void *p) {
      typedef vector<TGraph*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TGraph*>

namespace {
  void TriggerDictionaryInitialization_libSNsim_dict_Impl() {
    static const char* headers[] = {
"include/SNchannelNuE.hh",
"include/SNnupLS.hh",
"include/SNsource.hh",
"include/SNnumNakazatoSrc.hh",
"include/SNdetect.hh",
"include/SNchannelNCC.hh",
"include/SNpreGuoIntegFcn.hh",
"include/SNnumBurrowsSrc.hh",
"include/SNeffectLS.hh",
"include/SNnumGarchingSrc.hh",
"include/SNGarchingIntegFcn.hh",
"include/SNcncLS.hh",
"include/SNnccLS.hh",
"include/SNnumJapanSrc.hh",
"include/SNnumPreGuoSrc.hh",
"include/SNanaSrc.hh",
"include/SNchannelBCC.hh",
"include/SNchannelCNC.hh",
"include/SNnueLS.hh",
"include/SNpreKatoIntegFcn.hh",
"include/SNNakazatoIntegFcn.hh",
"include/SNchannelIBD.hh",
"include/SNchannelNuP.hh",
"include/SNchannels.hh",
"include/SNbccLS.hh",
"include/SNnumPreKatoSrc.hh",
"include/SNJapanIntegFcn.hh",
"include/SNibdLS.hh",
"include/SNBurrowsIntegFcn.hh",
0
    };
    static const char* includePaths[] = {
"/junofs/users/miaoyu/supernova/wenlj/simulation/include/",
"/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J21v1r0-Pre2/ExternalLibs/ROOT/6.22.08/include",
"/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J21v1r0-Pre1/ExternalLibs/ROOT/6.22.08/include/",
"/junofs/users/miaoyu/supernova/wenlj/simulation/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libSNsim_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$SNeffectLS.hh")))  __attribute__((annotate("$clingAutoload$include/SNchannelNuE.hh")))  SNeffectLS;
class __attribute__((annotate("$clingAutoload$SNchannels.hh")))  __attribute__((annotate("$clingAutoload$include/SNchannelNuE.hh")))  SNchannels;
class __attribute__((annotate("$clingAutoload$SNnueLS.hh")))  __attribute__((annotate("$clingAutoload$include/SNchannelNuE.hh")))  SNnueLS;
class __attribute__((annotate("$clingAutoload$include/SNchannelNuE.hh")))  SNchannelNuE;
class __attribute__((annotate("$clingAutoload$include/SNnupLS.hh")))  SNnupLS;
class __attribute__((annotate("$clingAutoload$include/SNsource.hh")))  SNsource;
class __attribute__((annotate("$clingAutoload$include/SNdetect.hh")))  SNdetect;
class __attribute__((annotate("$clingAutoload$SNnccLS.hh")))  __attribute__((annotate("$clingAutoload$include/SNchannelNCC.hh")))  SNnccLS;
class __attribute__((annotate("$clingAutoload$include/SNchannelNCC.hh")))  SNchannelNCC;
class __attribute__((annotate("$clingAutoload$include/SNpreGuoIntegFcn.hh")))  SNpreGuoIntegFcn;
class __attribute__((annotate("$clingAutoload$include/SNnumBurrowsSrc.hh")))  SNnumBurrowsSrc;
class __attribute__((annotate("$clingAutoload$include/SNnumGarchingSrc.hh")))  SNnumGarchingSrc;
class __attribute__((annotate("$clingAutoload$include/SNGarchingIntegFcn.hh")))  SNGarchingIntegFcn;
class __attribute__((annotate("$clingAutoload$include/SNcncLS.hh")))  SNcncLS;
class __attribute__((annotate("$clingAutoload$include/SNnumJapanSrc.hh")))  SNnumJapanSrc;
class __attribute__((annotate("$clingAutoload$include/SNnumPreGuoSrc.hh")))  SNnumPreGuoSrc;
class __attribute__((annotate("$clingAutoload$include/SNanaSrc.hh")))  SNanaSrc;
class __attribute__((annotate("$clingAutoload$SNbccLS.hh")))  __attribute__((annotate("$clingAutoload$include/SNchannelBCC.hh")))  SNbccLS;
class __attribute__((annotate("$clingAutoload$include/SNchannelBCC.hh")))  SNchannelBCC;
class __attribute__((annotate("$clingAutoload$include/SNchannelCNC.hh")))  SNchannelCNC;
class __attribute__((annotate("$clingAutoload$include/SNpreKatoIntegFcn.hh")))  SNpreKatoIntegFcn;
class __attribute__((annotate("$clingAutoload$SNibdLS.hh")))  __attribute__((annotate("$clingAutoload$include/SNchannelIBD.hh")))  SNibdLS;
class __attribute__((annotate("$clingAutoload$include/SNchannelIBD.hh")))  SNchannelIBD;
class __attribute__((annotate("$clingAutoload$include/SNchannelNuP.hh")))  SNchannelNuP;
class __attribute__((annotate("$clingAutoload$include/SNnumPreKatoSrc.hh")))  SNnumPreKatoSrc;
class __attribute__((annotate("$clingAutoload$include/SNJapanIntegFcn.hh")))  SNJapanIntegFcn;
class __attribute__((annotate("$clingAutoload$include/SNBurrowsIntegFcn.hh")))  SNBurrowsIntegFcn;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libSNsim_dict dictionary payload"

#ifndef _PGTRACK_
  #define _PGTRACK_ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "include/SNchannelNuE.hh"
#include "include/SNnupLS.hh"
#include "include/SNsource.hh"
#include "include/SNnumNakazatoSrc.hh"
#include "include/SNdetect.hh"
#include "include/SNchannelNCC.hh"
#include "include/SNpreGuoIntegFcn.hh"
#include "include/SNnumBurrowsSrc.hh"
#include "include/SNeffectLS.hh"
#include "include/SNnumGarchingSrc.hh"
#include "include/SNGarchingIntegFcn.hh"
#include "include/SNcncLS.hh"
#include "include/SNnccLS.hh"
#include "include/SNnumJapanSrc.hh"
#include "include/SNnumPreGuoSrc.hh"
#include "include/SNanaSrc.hh"
#include "include/SNchannelBCC.hh"
#include "include/SNchannelCNC.hh"
#include "include/SNnueLS.hh"
#include "include/SNpreKatoIntegFcn.hh"
#include "include/SNNakazatoIntegFcn.hh"
#include "include/SNchannelIBD.hh"
#include "include/SNchannelNuP.hh"
#include "include/SNchannels.hh"
#include "include/SNbccLS.hh"
#include "include/SNnumPreKatoSrc.hh"
#include "include/SNJapanIntegFcn.hh"
#include "include/SNibdLS.hh"
#include "include/SNBurrowsIntegFcn.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"SNBurrowsIntegFcn", payloadCode, "@",
"SNGarchingIntegFcn", payloadCode, "@",
"SNJapanIntegFcn", payloadCode, "@",
"SNanaSrc", payloadCode, "@",
"SNbccLS", payloadCode, "@",
"SNchannelBCC", payloadCode, "@",
"SNchannelCNC", payloadCode, "@",
"SNchannelIBD", payloadCode, "@",
"SNchannelNCC", payloadCode, "@",
"SNchannelNuE", payloadCode, "@",
"SNchannelNuP", payloadCode, "@",
"SNchannels", payloadCode, "@",
"SNcncLS", payloadCode, "@",
"SNdetect", payloadCode, "@",
"SNeffectLS", payloadCode, "@",
"SNibdLS", payloadCode, "@",
"SNnccLS", payloadCode, "@",
"SNnueLS", payloadCode, "@",
"SNnumBurrowsSrc", payloadCode, "@",
"SNnumGarchingSrc", payloadCode, "@",
"SNnumJapanSrc", payloadCode, "@",
"SNnumPreGuoSrc", payloadCode, "@",
"SNnumPreKatoSrc", payloadCode, "@",
"SNnupLS", payloadCode, "@",
"SNpreGuoIntegFcn", payloadCode, "@",
"SNpreKatoIntegFcn", payloadCode, "@",
"SNsource", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libSNsim_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libSNsim_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libSNsim_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libSNsim_dict() {
  TriggerDictionaryInitialization_libSNsim_dict_Impl();
}
