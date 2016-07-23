Blendenpik - A library for fast solution of rectangular systems.
Version 1.3, December 2009
Copyright (C) 2009, Haim Avron and Sivan Toledo.

Install
=======
Either FFTW or SPIRAL WHT are required. 

In MATLAB go to the directory of the unpacked files.
Type "install_blendenpik" which will guide you through the installation. 

Usage example
=============
>> A = rand(40000, 1000);
>> b = rand(40000, 1);
>> x = blendenpik(A, b);
		Random unit diagonal + unitary transformation time: 3.06 sec
		Random sampling time: 0.14 sec
		QR on random sample time: 2.28 sec
		Condition estimation: 0.44 sec
	Building preconditioner time: 6.00 sec
dense_lsqr: converged at iteration 42
	LSQR time: 7.20 sec
Total time: 13.22 sec

Paper
=====
The work is described in the following paper:

Haim Avron, Petar Maymounkov, and Sivan Toledo.
Blendenpik: Supercharging LAPACK's least-squares solver.
Submitted to SIAM Journal on Scientific Computing, August 2009.

Please cite it if you are using this software.

Copyright and License - BSD license
=====================
Copyright (c) 2009, Haim Avron and Sivan Toledo.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Tel-Aviv University nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY HAIM AVRON AND SIVAN TOLEDO ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

