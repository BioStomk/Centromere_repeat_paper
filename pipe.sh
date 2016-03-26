#!/bin/bash

SRCDIR=`pwd`/scripts
SPECIES=oryzias_latipes

export PERL5LIB=$SRCDIR/Perl_utils:$PERL5LIB

if [ ! -e logs ] ; then
	mkdir logs
fi


$SRCDIR/ftp_trace_reads.pl -one_species $SPECIES 2>&1 | tee logs/1_ftp_trace_reads.pl
$SRCDIR/parse_ftp_tracedb_data.pl -dir $SPECIES/ 2>&1 | tee logs/2_parse_ftp_tracedb_data.pl
$SRCDIR/trf_wrapper.pl -allfiles 2>&1 | tee logs/3_trf_wrapper.pl
$SRCDIR/trf_cluster.pl -trf $SPECIES.high.trf  2>&1 | tee logs/4_trf_cluster.pl
