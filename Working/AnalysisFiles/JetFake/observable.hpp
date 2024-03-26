//
//  observable.hpp
//  za
//
//  Created by Yandong Liu on 12/22/15.
//  Copyright © 2015 Yandong Liu. All rights reserved.
//

#ifndef observable_hpp
#define observable_hpp

#include <stdio.h>
#include <iostream>//needed for io
// #include "fastjet/ClusterSequence.hh"
using namespace fastjet;

struct sort_by_pt :
public std::binary_function<const PseudoJet&, const PseudoJet&, bool>
{
    bool operator()(const PseudoJet& x, const PseudoJet& y)
    {
        return  x.perp() > y.perp();
    }
};

struct sort_by_m :
public std::binary_function<const PseudoJet&, const PseudoJet&, bool>
{
    bool operator()(const PseudoJet& x, const PseudoJet& y)
    {
        return  x.m() > y.m();
    }
};

#endif /* observable_hpp */
