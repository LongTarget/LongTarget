#!/usr/bin/perl

use strict;
use Getopt::Std;
use Time::Local;
use DBI;
use Bio::EnsEMBL::Registry;

getopts('h:u:p:d:f:o:s:c:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $gtf_file, $output_path, $datasource_id, $species_id);

if ( defined($opts{h}) )
{
  $dbhost=$opts{h};
}
else
{
  $dbhost="localhost";
}

if ( defined($opts{u}) )
{
  $dbuser=$opts{u};
}
else
{
  $dbuser="master";
}

if ( defined($opts{o}) )
{
  $output_path=$opts{o};
}
else
{
  $output_path='./';
}

if ( defined($opts{p}))
{
  $dbpasswd=$opts{p};
}
else
{
  $dbpasswd='';
}

if ( defined($opts{d}))
{
  $dbname=$opts{d};
}
else
{
  $dbname='lncrna';
}

if ( defined($opts{f}) )
{
  $gtf_file=$opts{f};
}
else
{
  &usage();
  exit 1;
}

if ( defined($opts{s}))
{
  $datasource_id=$opts{s};
}
else
{
  $datasource_id=1;
}

if ( defined($opts{c}))
{
  $species_id=$opts{c};
}
else
{
  $species_id=1;
}

if (!open(GTFFILE, $gtf_file))
{
  die "Can't open $gtf_file.\n";
}

my $dbh = DBI->connect( "dbi:mysql:dbname=$dbname;host=$dbhost", $dbuser, $dbpasswd);
my $sth;
if ( !defined $dbh )
{
  close GTFFILE;
  die "Cannot connect to database  $dbhost $dbname $dbuser.\n";
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

if ($sth=$dbh->prepare("delete from tb_exon where datasource_id=$datasource_id"))
{
  $sth->execute();  
}
else
{
  $dbh->disconnect;
  close GTFFILE;
  die "Can't prepare SQL statement: $DBI::errstr\n";
}

if ($sth=$dbh->prepare("delete from tb_transcript where datasource_id=$datasource_id"))
{
  $sth->execute();
}
else
{
  $dbh->disconnect;
  close GTFFILE;
  die "Can't prepare SQL statement: $DBI::errstr\n";
}

if ($sth=$dbh->prepare("delete from tb_gene where datasource_id=$datasource_id"))
{
  $sth->execute();
}
else
{
  $dbh->disconnect;
  close GTFFILE;
  die "Can't prepare SQL statement: $DBI::errstr\n";
}

# get a slice adaptor for the human core database
my $species;
if( $species_id==1 )
{
  $species='Homo_sapiens';
}
else
{
  $species='Homo_sapiens';
}

my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );
my $slice;

my ($gene_id, $transcript_id);
my $gene_count=0;
my $transcript_count=0;
my $exon_count=0;
my $db_count=0;
while (<GTFFILE>)
{
  if ( /^#/ || /^\s+/ )
  {
    next;
  }

  chomp;
  my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attrs) = split /\t/;

  my $sqlstr;

  my ($gtf_gene_id,$gene_name,$havana_gene,$gene_status,$gene_type,$level,$sequence);
  my ($gtf_transcript_id,$transcript_name,$havana_transcript,$transcript_status,$transcript_type,$tag);
  my ($gtf_exon_id,$exon_number);

  if( $feature eq 'gene' )
  {
    $gtf_gene_id=&getAttr($attrs,'gene_id'); 
    $gene_name=&getAttr($attrs,'gene_name');
    $havana_gene=&getAttr($attrs,'havana_gene');
    $level=&getAttr($attrs,'level');
    $gene_status=&getAttr($attrs,'gene_status');
    $gene_type=&getAttr($attrs,'gene_type');
    my $chrom_num=substr($seqname,index($seqname,'chr')+3);
    $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom_num, $start, $end);
    if (!defined($slice)) 
    {
      $sequence='""';
    }
    else
    {
      $sequence = $slice->seq();
    }

    $sqlstr= "insert tb_gene (seqname,source,start,end,strand,gtf_gene_id,gene_name,havana_gene,gene_status,gene_type,level,datasource_id,species_id,sequence) value (\"$seqname\",\"$source\",$start,$end,\'$strand\',$gtf_gene_id,$gene_name,$havana_gene,$gene_status,$gene_type,$level,$datasource_id,$species_id,\'$sequence\')";
    #print "$sqlstr\n";
    $gene_count++;
  }
  elsif ( $feature eq 'transcript' )
  {
    $gtf_transcript_id=&getAttr($attrs,'transcript_id');
    $transcript_name=&getAttr($attrs,'transcript_name');
    $havana_transcript=&getAttr($attrs,'havana_transcript');
    $level=&getAttr($attrs,'level');
    $transcript_status=&getAttr($attrs,'transcript_status');
    $transcript_type=&getAttr($attrs,'transcript_type');
    $tag=&getAttr($attrs,'tag');
    $sqlstr= "insert tb_transcript (seqname,source,start,end,strand,gtf_transcript_id,transcript_name,havana_transcript,transcript_status,transcript_type,level,tag,datasource_id,gene_id) value (\"$seqname\",\"$source\",$start,$end,\'$strand\',$gtf_transcript_id,$transcript_name,$havana_transcript,$transcript_status,$transcript_type,$level,$tag,$datasource_id,$gene_id)";
    $transcript_count++;
  }
  elsif ( $feature eq 'exon' )
  {
    $gtf_exon_id=&getAttr($attrs,'exon_id');
    $exon_number=&getAttr($attrs,'exon_number');
    $level=&getAttr($attrs,'level');
    $transcript_status=&getAttr($attrs,'transcript_status');
    $transcript_type=&getAttr($attrs,'transcript_type');
    $tag=&getAttr($attrs,'tag');
    $sqlstr= "insert tb_exon (seqname,source,start,end,strand,gtf_exon_id,exon_number,transcript_status,transcript_type,level,tag,datasource_id,transcript_id) value (\"$seqname\",\"$source\",$start,$end,\'$strand\',$gtf_exon_id,$exon_number,$transcript_status,$transcript_type,$level,$tag,$datasource_id,$transcript_id)";
    $exon_count++;
  }
  else
  {
    print "Error in:$_\n";
    next;
  }
  if ($sth = $dbh->prepare( $sqlstr ))
  {
    #print "$sqlstr\n";
    if ($sth->execute())
    {
      $db_count++;
    }
    else
    {
      print "Error in:$_\n";
      print "Can't exec SQL statement: $DBI::errstr\n";
      next;
    }
  }
  else
  { 
    print "Can't prepare SQL statement: $DBI::errstr\n";
    next;
  }

  if ($sth = $dbh->prepare( 'select LAST_INSERT_ID()' ))
  {
    $sth->execute() or next;
    my ($last_id)=$sth->fetchrow_array();
    if( $feature eq 'gene' )
    {
      $gene_id=$last_id;
    }
    elsif ( $feature eq 'transcript' )
    {
      $transcript_id=$last_id;
    }
  }
  else
  {
    print "Can't prepare SQL statement: $DBI::errstr\n";
  }  
}

$sth->finish;

$dbh->disconnect;

close(GTFFILE);

print "Inserted gene records: $gene_count\n";
print "Inserted transcript records: $transcript_count\n";
print "Inserted exon records: $exon_count\n";
print "Inserted db records: $db_count\n";

sub getAttr
{
  my @fields=split /;\s*/, $_[0];
  my ($rec, @result);
  foreach $rec (@fields)
  {
    my $found=index($rec, $_[1]);
    if ( $found >= 0 )
    {
      @result=split /\s+/, $rec;
      #print "$result[0], $result[1]\n";
      return $result[1];
    }
  }
  if( $_[1] eq 'level' || $_[1] eq 'exon_number' )
  {
    return 0;
  }
  else
  {
    return '" "';
  }
}

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
  print "-f [GTF file name], default is '0'.\n";
  print "-o [Output result path].\n";
}

