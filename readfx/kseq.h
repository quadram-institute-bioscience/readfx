typedef void *readfx_gzfast_file;
int readfx_gzfast_read(readfx_gzfast_file fp, void *buf, int len);
#include "klib/kseq.h"
KSEQ_INIT(readfx_gzfast_file, readfx_gzfast_read)
