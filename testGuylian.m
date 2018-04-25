% addpath(genpath(pwd))

bitStream = CreateBitStream(768,1);
[bitsEnc, newH] = encoderLDPC_G(bitStream,128,256);

bitsRec = hardDecoderLPDC_G( bitsEnc, newH, 128, 256 );



% rmpath(genpath(pwd))



















