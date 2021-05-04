/*



--------- NUGA v1.0



*/
//Authors : S�m Landier (sam.landier@onera.fr)

#ifndef NUGA_HISTORY_H
#define NUGA_HISTORY_H

#include<vector>
#include "Nuga/include/defs.h"

namespace NUGA
{
  struct morse_t
  {
    std::vector<int> data, xr;
    
    morse_t() { xr.push_back(0); };
    
    void clear() { data.clear(); xr.clear(); xr.push_back(0); }
    void append(int val) { data.push_back(val); xr.push_back(data.size());/*next pos*/}
    void append(std::vector<int>& vals) { data.insert(data.end(), ALL(vals)); xr.push_back(data.size());}
  };

  struct history_t
  {
    std::vector<int> ptoids, ptnids;
    morse_t pgoids, phoids, pgnids, phnids;

    void clear() { ptoids.clear();  pgoids.clear();  phoids.clear();  ptnids.clear(); pgnids.clear();  phnids.clear(); }

    template<typename T>
    void transfer_pg_colors(std::vector<E_Int>& src_ids, std::vector<T> & src_flags, std::vector<T>& target_flags)
    {
      //compute target size : find max id
      E_Int target_sz = 0;
      for (size_t i=0; i < src_ids.size(); ++i)
      {
        E_Int srcid = src_ids[i];

        E_Int nb_target = pgnids.xr[srcid+1] - pgnids.xr[srcid];
        const int* tgt_start = &pgnids.data[pgnids.xr[srcid]];
        for (size_t t=0; t < nb_target; ++t)
        {
          if (tgt_start[t] == IDX_NONE) continue;
          if (tgt_start[t] < 0) target_sz = std::max(target_sz, -(tgt_start[t]+1));
          else target_sz = std::max(target_sz, tgt_start[t]);
        }
      }

      target_flags.resize(target_sz, IDX_NONE);

      for (size_t i=0; i < src_ids.size(); ++i)
      {
        E_Int srcid = src_ids[i];

        E_Int nb_target = pgnids.xr[srcid+1] - pgnids.xr[srcid];
        const int* tgt_start = &pgnids.data[pgnids.xr[srcid]];

        if (nb_target == 1)
        {
          if (*tgt_start == IDX_NONE) // is gone
            continue;
          else if (*tgt_start < 0)    //  is agglomerated into target
          {
            E_Int tgtid = -(*tgt_start+1);
            if (target_flags[tgtid] == -1) //frozen because contrinutions have different flags
              continue;
            else if (target_flags[tgtid] == IDX_NONE)
              target_flags[tgtid] = src_flags[srcid];
            else if (target_flags[tgtid] != src_flags[srcid]) //freeze
              target_flags[tgtid] = -1;
          }
          else // just a move 
            target_flags[*tgt_start] = src_flags[srcid];
        }
        else // it's a split : all children inherits parent value
        {
          for (size_t t=0; t < nb_target; ++t)
            target_flags[tgt_start[t]] = src_flags[srcid];
        }
      }
    }

  }; //history_t
}

#endif