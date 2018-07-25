#!/usr/bin/perl

use strict;
use Getopt::Std;
use Time::Local;
use DBI;

getopts('h:u:p:d:f:s:c:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $gtf_file, $datasource_id, $species_id);

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

if ($sth=$dbh->prepare("delete from tb_myexon where datasource_id=$datasource_id"))
{
  $sth->execute();  
}
else
{
  $dbh->disconnect;
  close GTFFILE;
  die "Can't prepare SQL statement: $DBI::errstr\n";
}

my ($gene_id,$db_count,$exon_number);
my $db_count=0;
while (<GTFFILE>)
{
  if ( /^#/ || /^\s+/ )
  {
    next;
  }

  chomp;
  my ($seqname, $feature, $start, $end, $strand, $gtf_gene_id, $gtf_exon_id) = split /\t/;

  my $sqlstr;

  if( $feature eq 'gene' )
  {
    $exon_number=0;
    $sqlstr= "select gene_id from tb_gene where gtf_gene_id=\'$gtf_gene_id\'";
    if ($sth = $dbh->prepare( $sqlstr ))
    {
      if ($sth->execute())
      {
        if(!(($gene_id)=$sth->fetchrow_array()))
        {
          $gene_id=0;
        } 
      }
      else
      {
        print "Error in:$_\n";
        print "Can't exec SQL statement: $DBI::errstr\n";
        $gene_id=0;next;
      }
    }
    else
    {
      print "Can't prepare SQL statement: $DBI::errstr\n";
      next;
    }
  }
  elsif ( ($feature =~ /^exon/)  && ($gene_id > 0) )
  {
    $exon_number++;
    $sqlstr= "insert tb_myexon (seqname,start,end,strand,gtf_exon_id,exon_number,datasource_id,gene_id) value (\"$seqname\",$start,$end,\'$strand\',\'$gtf_exon_id\',$exon_number,$datasource_id,$gene_id)";
    if ($sth = $dbh->prepare( $sqlstr ))
    {
      if ($sth->execute())
      {
        $db_count++; 
      }
      else
      {
        print "Error in:$_\n";
        print "Can't exec SQL statement: $DBI::errstr\n";
        $gene_id=0;next;
      }
    }
    else
    {
      print "Can't prepare SQL statement: $DBI::errstr\n";
      next;
    }
  }
  else
  {
    print "Error in:$_\n";
    next;
  }
  print "Writen db records: $db_count\r";
}
print "\n";
$sth->finish;

$dbh->disconnect;

close(GTFFILE);

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
  print "-f [GTF file name], default is '0'.\n";
}

