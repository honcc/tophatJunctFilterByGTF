#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use List::Util qw ( sum );
use Storable;
use Data::Dumper::Names;

######################################################################################################################################################
#
#	Description
#		This is a perl script to extract the junctions in a tophat BED file based on the junctions in a GTF file. It reads th BED file junctions and junctions in the GTF 
#	, then it prints the junctions overlapping between the two.
#
#	Input
#		--BEDPathListFile=		file path; path of the a file contains the list of the input BED files;
#		--GTFPath=				file path; path of the GTF file;
#		--BEDFormat=			string; 'tophat' or 'HMMSplicerBEDToSAMParser'; To determine which column to get the readnum; default=tophat;
#		--cuffDiffGenePath=		file path; gene_exp.diff output from cuffdiff for read of the sample gene expression level and differential expression;
#		--cuffDiffIsoformPath=	file path; isoform_exp.diff output from cuffdiff for read of the sample gene expression level and differential expression;
#		--outDir=				dir path: path of the output dir;
#
#	Usage
#		
#		perl tophatJunctFilterByGTF_v0.1.pl --BEDPathListFile=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/splicing/tophatJunctFilterByGTF/v0.1/HM1Rahman.txt --GTFPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.bfmRNA.nonStochasticIsoform.gtf --outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/splicing/tophatJunctFilterByGTF/v0.1/tophatJunctFilterByGTF/
#
#	Assumption
#
#	History:
#		
#		v0.1
#		-debut
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#
#----------Read parameters ----------#
my ($BEDPathListFile, $GTFPath, $BEDFormat, $cuffDiffGenePath, $cuffDiffIsoformPath, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

my $geneCuffdiffDataHsh_ref = readCuffdiffData($cuffDiffGenePath);
my $isoformCuffdiffDataHsh_ref = readCuffdiffData($cuffDiffIsoformPath);

my ($BEDInfoHsh_ref, $cuffdiffSampleNumHsh_ref) = readBEDFilePath($BEDPathListFile);

my $BEDDataHsh_ref;
($BEDDataHsh_ref, $BEDInfoHsh_ref) = readAllBEDFile($BEDInfoHsh_ref, $BEDFormat);
my %BEDInfoHsh = %{$BEDInfoHsh_ref};

my ($GTFJunctHsh_ref, $GTFDataHsh_ref) = getGTFJunctStr($GTFPath);

my $cutoffSignificantRdNum = 5;
filterBEDJunction($BEDDataHsh_ref, $GTFJunctHsh_ref, \%BEDInfoHsh, $GTFDataHsh_ref, $BEDFormat, $outDir, $cutoffSignificantRdNum, $geneCuffdiffDataHsh_ref, $isoformCuffdiffDataHsh_ref, $cuffdiffSampleNumHsh_ref);

printCMDLogOrFinishMessage("finishMessage");

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	my ($BEDPathListFile, $GTFPath, $BEDFormat, $cuffDiffGenePath, $cuffDiffIsoformPath, $outDir);
	
	$outDir = "./tophatJunctFilterByGTF/";
	$BEDFormat = 'topHat';

	foreach my $param (@ARGV) {
		if ($param =~ m/--BEDPathListFile=/) {$BEDPathListFile = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--GTFPath=/) {$GTFPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--BEDFormat=/) {$BEDFormat = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--cuffDiffGenePath=/) {$cuffDiffGenePath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--cuffDiffIsoformPath=/) {$cuffDiffIsoformPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);} 
	}
	
	#---check file
	foreach my $fileToCheck ($BEDPathListFile, $GTFPath, $cuffDiffGenePath, $cuffDiffIsoformPath) {
		die "Cant read $fileToCheck" if not -s $fileToCheck;
		my @filePathAry = split /\//, $fileToCheck;
		my $fileName = $filePathAry[-1];
		print "$fileName checked.\n";
	}

	my @paramAry = ($BEDPathListFile, $GTFPath, $BEDFormat, $cuffDiffGenePath, $cuffDiffIsoformPath, $outDir);

	system "mkdir -p -m 777 $outDir/";
	open (PARAM, ">$outDir/parameters.txt");
	print PARAM Dumper($BEDPathListFile, $GTFPath, $BEDFormat, $cuffDiffGenePath, $cuffDiffIsoformPath, $outDir);
	close PARAM;
	
	return @paramAry;
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## readBEDFilePath
sub readBEDFilePath {

	#---my $BEDInfoHsh_ref = readBEDFilePath($BEDPathListFile);
	my $BEDPathListFile = $_[0];
	
	my %BEDInfoHsh;
	my %cuffdiffSampleNumHsh;
	open (BEDPATHFILE, "$BEDPathListFile");
	while (my $theLine = <BEDPATHFILE>) {
		chomp $theLine;
		my ($cuffdiffSampleNum, $sampleName, $replicate, $BEDPath) = split /\t/, $theLine;
		$cuffdiffSampleNumHsh{$cuffdiffSampleNum} = $sampleName;
		${${$BEDInfoHsh{$sampleName}}{$replicate}}{'BEDPath'} = $BEDPath;
		die "Cant read $BEDPath" if not -s $BEDPath;
		print "$sampleName $replicate BED checked.\n";
	}
	close BEDPATHFILE;
	
	return \%BEDInfoHsh, \%cuffdiffSampleNumHsh;
}
########################################################################## readIndividualBEDFile
sub readIndividualBEDFile {

	#---my $BEDInfoHsh_ref = readIndividualBEDFile($BEDPath);
	
	my $BEDPath = $_[0];
	my $BEDFormat = $_[1];
	
	my %indivBEDDataHsh;
	my $junctNum = 0;
	my $totalReadNum = 0;
	open (BED, $BEDPath);
	while (my $theLine = <BED>) {
		chomp $theLine;
		next if ($theLine =~ m/^track name/);
		$junctNum++;
		my @theLineSplt = split (/\t/, $theLine);
		#DS571145	13916	14145	JUNC00000008	141	-	13916	14145	255,0,0	2	77,97	0,132
		my $cntg = $theLineSplt[0];
		my $bedStart = $theLineSplt[1];
		my $bedEnd = $theLineSplt[2];
		my $readNum;
		if ($BEDFormat eq 'HMMSplicerBEDToSAMParser') {
			$readNum = $1 if $theLineSplt[3] =~ m/read=(\d+)\|/;
			die 'BED format is not HMMSplicerBEDToSAMParser. Quitting' if not defined $readNum;
		} elsif ($BEDFormat eq 'topHat') {
			$readNum = $theLineSplt[4];
		}
		$totalReadNum += $readNum;
		my $strd = $theLineSplt[5]; #--- will be an asterick for non-cannoical
		my @blkSizesSplt = split /,/, $theLineSplt[10];
		my $blk1Size = $blkSizesSplt[0];
		my $blk2Size = $blkSizesSplt[1];
		
		my $intronStart = $bedStart + $blk1Size + 1;
		my $intronEnd = $bedEnd - $blk2Size;

		my $junctStr = $cntg.":".$intronStart.":".$intronEnd; #---assumed to be unique
		
		${$indivBEDDataHsh{$junctStr}}{'strnd'} = $strd;
		${$indivBEDDataHsh{$junctStr}}{'readNum'} = $readNum;
		my $BEDLine = $theLine;
		${$indivBEDDataHsh{$junctStr}}{'BEDLine'} = $BEDLine;

	}
	close BED;
	
	print "Totally $junctNum junctions read.\n";
	
	return \%indivBEDDataHsh, $totalReadNum;
}
########################################################################## readAllBEDFile
sub readAllBEDFile {
	
	#---my $BEDDataHsh_ref = readAllBEDFile($BEDInfoHsh_ref);

	my %BEDInfoHsh = %{$_[0]};
	my $BEDFormat = $_[1];
	
	my %BEDDataHsh;
	my @allTotalReadNumAry;
	foreach my $sampleName (sort {$a cmp $b} keys %BEDInfoHsh) {
		foreach my $replicate (sort {$a <=> $b} keys %{$BEDInfoHsh{$sampleName}}) {
			print "Reading BED of $sampleName $replicate.\n";
			my $BEDPath = ${${$BEDInfoHsh{$sampleName}}{$replicate}}{'BEDPath'};
			my ($indivBEDDataHsh_ref, $totalReadNum) = readIndividualBEDFile($BEDPath, $BEDFormat);
			${$BEDInfoHsh{$sampleName}{$replicate}}{'totalReadNum'} = $totalReadNum;
			push @allTotalReadNumAry, $totalReadNum;
			%{${$BEDDataHsh{$sampleName}}{$replicate}} = %{$indivBEDDataHsh_ref};
		}
	}
	
	my $avgTotalRdNum = sum(@allTotalReadNumAry)/@allTotalReadNumAry;
	
	foreach my $sampleName (sort {$a cmp $b} keys %BEDInfoHsh) {
		foreach my $replicate (sort {$a <=> $b} keys %{$BEDInfoHsh{$sampleName}}) {
			my $totalReadNum = ${$BEDInfoHsh{$sampleName}{$replicate}}{'totalReadNum'};
			my $scaleFactor = $avgTotalRdNum/$totalReadNum;
			${$BEDInfoHsh{$sampleName}{$replicate}}{'scaleFactor'} = $scaleFactor;
			print "$sampleName $replicate: totalReadNum = $totalReadNum, scaleFactor = $scaleFactor\n"
		}
	}

	return \%BEDDataHsh, \%BEDInfoHsh;
}
########################################################################## getGTFJunctStr
sub getGTFJunctStr {
	
	#---my ($GTFJunctHsh_ref, $GTFDataHsh_ref) = getGTFJunctStr($GTFPath);
	
	my $GTFPath = $_[0];
	
	my %GTFDataHsh;
	my %transcriptCountHsh;
	my %GTFJunctHsh;
	
	my %geneJunctStrHsh;

	open (GTF, "$GTFPath");
	while (my $theLine = <GTF>) {
		chomp $theLine;
		#DS571190	IPasteur	exon	7775	8063	.	+	.	gene_id "EHI_038570"; transcript_id "EHI_038570.ref"; gene_name "RNA polymerase subunit Rpb8"; transcript_name "RNA polymerase subunit Rpb8";
		my ($sampleName, $BEDPath) = split /\t/, $theLine;
		my ($cntg, $source, $type, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $theLine;
		my $gene_id = $1 if $attribute =~ m/gene_id "(\S+)";/;
		my $transcript_id = $1 if $attribute =~ m/transcript_id "(\S+)";/;
		my $gene_name = $1 if $attribute =~ m/gene_name "(\S+)";/;
		my $transcript_name = $1 if $attribute =~ m/transcript_name "(\S+)";/;
		
		push @{${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'boundary'}}, ($start, $end);
		${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'strand'} = $strand;
		${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'transcript_name'} = $transcript_name;
		${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'cntg'} = $cntg;
		
		$transcriptCountHsh{$transcript_id}++;
	}
	
	my $geneNum = keys %GTFDataHsh;
	my $transciptNum = keys %transcriptCountHsh;
	
	print "$geneNum gene and $transciptNum transcripts read.\n";
	
	foreach my $gene_id (keys %GTFDataHsh) {
		foreach my $transcript_id (keys %{$GTFDataHsh{$gene_id}}) {
			next if @{${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'boundary'}} == 2; #---no junction
			my @sortedBoundAry = sort {$a <=> $b} @{${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'boundary'}};
			my $cntg = ${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'cntg'};
			for (my $i = 1; $i < $#sortedBoundAry; $i=$i+2) {
				my $intronStart = $sortedBoundAry[$i]+1;
				my $intronEnd = $sortedBoundAry[$i+1]-1;
				my $junctStr = $cntg.":".$intronStart.":".$intronEnd; #---assumed to be unique
				${$GTFJunctHsh{$junctStr}}{'gene_id'} = $gene_id;
				push @{${$GTFJunctHsh{$junctStr}}{'transcript_id'}}, $transcript_id;
				push @{${${$GTFDataHsh{$gene_id}}{$transcript_id}}{'junctStr'}}, $junctStr;

			}
		}
	}
	my $junctNum = keys %GTFJunctHsh;
	
	print "Totally $junctNum junctions stored from GTF.\n";
	
	return \%GTFJunctHsh, \%GTFDataHsh;

}
########################################################################## filterBEDJunction
sub filterBEDJunction {
	
	#---filterBEDJunction($BEDDataHsh_ref, $GTFJunctHsh_ref);
	my %BEDDataHsh = %{$_[0]};
	my %GTFJunctHsh = %{$_[1]};
	my %BEDInfoHsh = %{$_[2]};
	my %GTFDataHsh = %{$_[3]};
	my $BEDFormat = $_[4];
	my $outDir = $_[5];
	my $cutoffSignificantRdNum = $_[6];
	my %geneCuffdiffDataHsh = %{$_[7]};
	my %isoformCuffdiffDataHsh = %{$_[8]};
	my %cuffdiffSampleNumHsh = %{$_[9]};
	
	my $blk1Size = my $blk2Size = 5;
	my %junctCountHsh;
	
	foreach my $sampleName (sort {$a cmp $b} keys %BEDDataHsh) {
		foreach my $replicate (sort {$a <=> $b} keys %{$BEDDataHsh{$sampleName}}) {
			my $junctHit = 0;
			my $scaleFactor = ${$BEDInfoHsh{$sampleName}{$replicate}}{'scaleFactor'};
			open (ORIGINALFILTER, ">$outDir/arc.original.filter.$sampleName.$replicate.BED");
			open (SCALEDFILTERARC, ">$outDir/arc.scaled.filter.$sampleName.$replicate.BED");
			open (SCALEDFILTERBOX, ">$outDir/box.scaled.filter.$sampleName.$replicate.BED");
			print ORIGINALFILTER "track name=junctions description=\"extracted TopHat junctions\"\n";
			print SCALEDFILTERARC "track name=junctions description=\"extracted TopHat junctions\"\n";
			print SCALEDFILTERBOX "track name=$sampleName.$replicate description=\"extracted TopHat junctions\"\n";
	
			foreach my $junctStr (sort {$a cmp $b} keys %{$BEDDataHsh{$sampleName}{$replicate}}) {
				
				next if not exists $GTFJunctHsh{$junctStr};
				$junctHit++;
				my $BEDLine = ${${$BEDDataHsh{$sampleName}{$replicate}}{$junctStr}}{'BEDLine'};
				my @BEDLineSplt = split /\t/, $BEDLine;

				my $readNum;
				if ($BEDFormat eq 'HMMSplicerBEDToSAMParser') {
					$readNum = $1 if $BEDLineSplt[3] =~ m/read=(\d+)\|/;
					die 'BED format is not HMMSplicerBEDToSAMParser. Quitting' if not defined $readNum;
				} elsif ($BEDFormat eq 'topHat') {
					$readNum = $BEDLineSplt[4];
				}
	
				my $scaledReadNum = sprintf "%.2f", ($readNum*$scaleFactor);
				$BEDLineSplt[3] = "read(scaled)=$readNum($scaledReadNum)";
				$BEDLineSplt[4] = sprintf "%.0f", $scaledReadNum;
				my ($cntg, $intronStart, $intronEnd) = split /:/, $junctStr;
				my $bedStart = $intronStart - $blk1Size - 1;
				my $bedEnd = $intronEnd + $blk2Size;
				$BEDLineSplt[1] = $bedStart;
				$BEDLineSplt[2] = $bedEnd;
				my $blk1Start = 0;
				my $blk2Start = $blk1Size + ($intronEnd - $intronStart) + 1;
				$BEDLineSplt[10] = "$blk1Size,$blk2Size";
				$BEDLineSplt[11] = "$blk1Start,$blk2Start";
	
				my $scaledBEDLine = join "\t", @BEDLineSplt;
				print ORIGINALFILTER $BEDLine."\n";
				print SCALEDFILTERARC $scaledBEDLine."\n";
				print SCALEDFILTERBOX $scaledBEDLine."\n";

				#---collect the junctCount data
				${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'readNum'} = $readNum;
				${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'scaledReadNum'} = $scaledReadNum;

			}
			
			print "$junctHit junctions were matched in $sampleName.\n";
			close ORIGINALFILTER;
			close SCALEDFILTERBOX;
			close SCALEDFILTERARC;
		}
	}
	
	#---fill zeros
	foreach my $junctStr (keys %GTFJunctHsh) { 
		foreach my $sampleName (sort {$a cmp $b} keys %BEDDataHsh) {
			foreach my $replicate (sort {$a <=> $b} keys %{$BEDDataHsh{$sampleName}}) {
				if (not exists ${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}) {
					${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'readNum'} = 0;
					${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'scaledReadNum'} = 0;
				}
			}
		}
	}

	my %nonRedundantSamplePairHsh;
	foreach my $sampleNameRef (sort {$a cmp $b} keys %BEDDataHsh) {
		foreach my $sampleNameQry (sort {$a cmp $b} keys %BEDDataHsh) {
			next if $sampleNameRef eq $sampleNameQry;
			next if exists ${$nonRedundantSamplePairHsh{$sampleNameRef}}{$sampleNameQry};
			next if exists ${$nonRedundantSamplePairHsh{$sampleNameQry}}{$sampleNameRef};
			${$nonRedundantSamplePairHsh{$sampleNameRef}}{$sampleNameQry}++;
		}
	}

	
	open (JUNCTSAMPLECOUNT, ">$outDir/junctionBasedCountBetweenSamples.log.xls");
	foreach my $junctStr (keys %junctCountHsh) { 
		my @outputAry;
		push @outputAry, 'junctStr';
		push @outputAry, 'gene_id';
		push @outputAry, 'geneName';
		push @outputAry, 'numTrnscrptInGene';
		push @outputAry, 'numTrnscrptForJunct';
		push @outputAry, 'transcript_id';
		push @outputAry, 'trnscrptSpecific';
		push @outputAry, 'refOrAlt';

		foreach my $cuffdiffSampleNum (sort {$a cmp $b} keys %cuffdiffSampleNumHsh) {
			my $sampleName = $cuffdiffSampleNumHsh{$cuffdiffSampleNum};
			push @outputAry, "gene.exp.$sampleName";
		}
		push @outputAry, "gene.exp.significant";

		foreach my $cuffdiffSampleNum (sort {$a cmp $b} keys %cuffdiffSampleNumHsh) {
			my $sampleName = $cuffdiffSampleNumHsh{$cuffdiffSampleNum};
			push @outputAry, "isfm.exp.$sampleName";
		}

		push @outputAry, "isfmpVal";
		push @outputAry, "isfmqVal";
		push @outputAry, "isfmLog2FC";
		push @outputAry, "isfmSignificant";
		
		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
				push @outputAry, "readNum.$sampleName.$replicate";
			}
		}
		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
			}
			push @outputAry, "readNum.$sampleName.Total";
		}

		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
				push @outputAry, "scaledReadNum.$sampleName.$replicate";
			}
		}
		
		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
			}
			push @outputAry, "scaledReadNum.$sampleName.Total";
		}
		
		foreach my $sampleNameRef (sort {$a cmp $b} keys %nonRedundantSamplePairHsh) {
			foreach my $sampleNameQry (sort {$a cmp $b} keys %{$nonRedundantSamplePairHsh{$sampleNameRef}}) {
				push @outputAry, "$sampleNameRef.vs.$sampleNameQry.scaledReadNum.ratio";
			}
		}

		push @outputAry, "allSampleRdNum>$cutoffSignificantRdNum";
		
		print JUNCTSAMPLECOUNT join '', ((join "\t", @outputAry), "\n");
		last;
	}
	
	my %sampleBasedCountHsh;
	
	foreach my $junctStr (keys %junctCountHsh) { 
		my @outputAry;
		push @outputAry, $junctStr;
		my $gene_id = ${$GTFJunctHsh{$junctStr}}{'gene_id'};
		my $numTrnscrptForJunct = @{${$GTFJunctHsh{$junctStr}}{'transcript_id'}};
		my $numTrnscrptInGene = keys %{$GTFDataHsh{$gene_id}};
		my $transcript_id = join ";", @{${$GTFJunctHsh{$junctStr}}{'transcript_id'}};
		my $trnscrptSpecific = 'no';
		$trnscrptSpecific = 'yes' if (($numTrnscrptForJunct < $numTrnscrptInGene) and ($numTrnscrptForJunct == 1));
		
		my $refOrAlt = 'alt';
		$refOrAlt = 'ref' if $transcript_id =~ m/\.ref$/;
		my $geneName = ${$geneCuffdiffDataHsh{$gene_id}}{'geneName'};
		
		push @outputAry, $gene_id;
		push @outputAry, $geneName;
		push @outputAry, $numTrnscrptInGene;
		push @outputAry, $numTrnscrptForJunct;
		push @outputAry, $transcript_id;
		push @outputAry, $trnscrptSpecific;
		push @outputAry, $refOrAlt;

		foreach my $cuffdiffSampleNum (sort {$a cmp $b} keys %cuffdiffSampleNumHsh) {
			my $geneRPKM = ${$geneCuffdiffDataHsh{$gene_id}}{$cuffdiffSampleNum};
			push @outputAry, $geneRPKM;
		}
		
		my $geneSignificant = ${$geneCuffdiffDataHsh{$gene_id}}{'significant'};
		push @outputAry, $geneSignificant;

		foreach my $cuffdiffSampleNum (sort {$a cmp $b} keys %cuffdiffSampleNumHsh) {
			my $sampleName = $cuffdiffSampleNumHsh{$cuffdiffSampleNum};
			my $isfmRPKM = 'null';
			$isfmRPKM = ${$isoformCuffdiffDataHsh{$transcript_id}}{$cuffdiffSampleNum} if $trnscrptSpecific eq 'yes';
			push @outputAry, $isfmRPKM;
		}
		
		my $isfmSignificant = my $isfmLog2FC = my $isfmpVal = my $isfmqVal = 'null';
		if ($trnscrptSpecific eq 'yes') {
			$isfmLog2FC = ${$isoformCuffdiffDataHsh{$transcript_id}}{'log2fold_change'};
			$isfmpVal = ${$isoformCuffdiffDataHsh{$transcript_id}}{'p_value'};
			$isfmqVal = ${$isoformCuffdiffDataHsh{$transcript_id}}{'q_value'};
			$isfmSignificant = ${$isoformCuffdiffDataHsh{$transcript_id}}{'significant'};
		}
		
		push @outputAry, $isfmpVal;
		push @outputAry, $isfmqVal;
		push @outputAry, $isfmLog2FC;
		push @outputAry, $isfmSignificant;

		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
				push @outputAry, ${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'readNum'};
			}
		}
		my $allSampleRdNumSignificant = 'yes';
		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			my $totalNum = 0;
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
				$totalNum += ${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'readNum'};
			}
			${${$sampleBasedCountHsh{$junctStr}}{'readNum'}}{$sampleName} = $totalNum;
			$allSampleRdNumSignificant = 'no' if $totalNum < $cutoffSignificantRdNum;
			push @outputAry, $totalNum;
		}

		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
				push @outputAry, ${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'scaledReadNum'};
			}
		}
		
		foreach my $sampleName (sort {$a cmp $b} keys %{$junctCountHsh{$junctStr}}) {
			my $totalNum = 0;
			foreach my $replicate (sort {$a <=> $b} keys %{${$junctCountHsh{$junctStr}}{$sampleName}}) {
				$totalNum += ${${${$junctCountHsh{$junctStr}}{$sampleName}}{$replicate}}{'scaledReadNum'};
			}
			${${$sampleBasedCountHsh{$junctStr}}{'scaledReadNum'}}{$sampleName} = $totalNum;
			push @outputAry, $totalNum;
		}
		
		#---ratio of scaledReadNum
		foreach my $sampleNameRef (sort {$a cmp $b} keys %nonRedundantSamplePairHsh) {
			foreach my $sampleNameQry (sort {$a cmp $b} keys %{$nonRedundantSamplePairHsh{$sampleNameRef}}) {
				my $refScaledReadNum = ${${$sampleBasedCountHsh{$junctStr}}{'scaledReadNum'}}{$sampleNameRef};
				my $qryScaledReadNum = ${${$sampleBasedCountHsh{$junctStr}}{'scaledReadNum'}}{$sampleNameQry};
				my $log2Ratio;
				if (($qryScaledReadNum > 0) and ($refScaledReadNum > 0)) {
					$log2Ratio = log($qryScaledReadNum/$refScaledReadNum)/log(2);

				} elsif (($qryScaledReadNum == 0) and ($refScaledReadNum == 0)) {
					$log2Ratio = 'null';
				
				} elsif (($qryScaledReadNum > 0) and ($refScaledReadNum == 0)) {
					$log2Ratio = 9999;
				
				} elsif (($qryScaledReadNum == 0) and ($refScaledReadNum > 0)) {
					$log2Ratio = -9999;

				} else {
				
				}
				push @outputAry, $log2Ratio;
			}
		}
		
		push @outputAry, $allSampleRdNumSignificant;
		
		print JUNCTSAMPLECOUNT join '', ((join "\t", @outputAry), "\n");
		
	}
	close JUNCTSAMPLECOUNT;
}
########################################################################## readCuffdiffData
sub readCuffdiffData {
	
	my $cuffDiffPath = $_[0];
	
	print "Reading $cuffDiffPath.\n";
	
	my %cuffDiffDataHsh;
	open (CUFFDIFF, "$cuffDiffPath");
	my $header = <CUFFDIFF>;
	while (my $theLine = <CUFFDIFF>) {
		chomp $theLine;
		my ($test_id, $gene_id, $gene, $locus, $sample_1, $sample_2, $status, $value_1, $value_2, $log2fold_change, $test_stat, $p_value, $q_value, $significant) = split /\t/, $theLine;
		${$cuffDiffDataHsh{$test_id}}{'gene_id'} = $gene_id;
		${$cuffDiffDataHsh{$test_id}}{'geneName'} = $gene;
		${$cuffDiffDataHsh{$test_id}}{'value_1'} = $value_1;
		${$cuffDiffDataHsh{$test_id}}{'value_2'} = $value_2;
		${$cuffDiffDataHsh{$test_id}}{'log2fold_change'} = $log2fold_change;
		${$cuffDiffDataHsh{$test_id}}{'q_value'} = $q_value;
		${$cuffDiffDataHsh{$test_id}}{'p_value'} = $p_value;
		${$cuffDiffDataHsh{$test_id}}{'significant'} = $significant;
	}
	close (CUFFDIFF);
	
	my $itemNum = keys %cuffDiffDataHsh;
	
	print "cuffdiff data of $itemNum items stored.\n";
	
	return \%cuffDiffDataHsh;
}
