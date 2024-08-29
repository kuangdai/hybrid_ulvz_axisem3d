//
//  ftz.c
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/6/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  To activate flush to zero for denormal float handling
//  http://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x

#ifndef _CRAYC
#include <xmmintrin.h>
#endif

void set_ftz() {
#ifndef _CRAYC
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#ifndef NDEBUG
    _MM_SET_EXCEPTION_MASK(_MM_MASK_INEXACT |
                           _MM_MASK_UNDERFLOW |
                           _MM_MASK_DENORM);
#endif
#endif
}
