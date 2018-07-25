#!/usr/bin/perl

use strict;
use Getopt::Std;
use DBI;

getopts('h:u:p:d:s:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $datasource_id);

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
  $dbname='lncrnadb';
}

if ( defined($opts{s}))
{
  $datasource_id=$opts{s};
}
else
{
  $datasource_id=1;
}

my $dbh = DBI->connect( "dbi:mysql:dbname=$dbname;host=$dbhost", $dbuser, $dbpasswd);
my ($sth1, $sth2, $sth3);
if ( !defined $dbh )
{
  die "Cannot connect to database  $dbhost $dbname $dbuser.\n";
}

if ($sth1=$dbh->prepare("delete from tb_splicer where datasource_id=$datasource_id"))
{
  $sth1->execute();
}
else
{
  $dbh->disconnect;
  die "Can't prepare SQL statement: $DBI::errstr\n";
}

my $gene_count=0;
my $exon_count=0;
my $exon_count_max=0;
my $db_count=0;
my $head=0;
my $tail=0;
my $head_splicer='';
my $tail_splicer='';

$sth1 = $dbh->prepare("select gene_id,start,end,flank,sequence from tb_gene where datasource_id=$datasource_id")
  or die "Can't prepare SQL statement: $DBI::errstr\n";
$sth1->execute()  or die "Can't execute SQL statement: $DBI::errstr\n";

my($gene_id, $start, $end, $flank, $sequence);
my($exon_id,$last_exon_id,$exon_start,$exon_end,$exon_strand,$exon_number);
while(($gene_id, $start, $end, $flank, $sequence) = $sth1->fetchrow_array())
{
  $gene_count++;

  if ($sth2 = $dbh->prepare( "select count(*) from tb_myexon where gene_id=$gene_id"))
  {
    if( $sth2->execute())
    {
      if( !(($exon_count_max) = $sth2->fetchrow_array()) )
      {
        print "\nCan't fetchrow_array: $DBI::errstr\n";
        next;
      }
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
    print "\nCan't prepare SQL statement: $DBI::errstr\n";
    next;
  } 

  if ($sth2 = $dbh->prepare( "select exon_id,start,end,strand,exon_number from tb_myexon where gene_id=$gene_id order by exon_number"))
  {
    if( $sth2->execute())
    {
      while( ($exon_id,$exon_start,$exon_end,$exon_strand,$exon_number) = $sth2->fetchrow_array() )
      {
        $exon_count++;
        $last_exon_id=$exon_id;

        #Check head splicer
        if ($exon_strand eq '+')
        {
          $head_splicer = substr($sequence, $exon_start+$flank-$start-2, 2);substr($sequence, $exon_start-$start-2, 2);
          $tail_splicer = substr($sequence, $exon_end+$flank-$start+1, 2);substr($sequence, $exon_end-$start, 2);
        }
        else
        {
          $head_splicer = substr($sequence, $end+$flank-$exon_end-2, 2);
          $tail_splicer = substr($sequence, $end+$flank-$exon_start+1, 2);
        }

        if( $exon_number == 1 )
        {
          $head=0;
        }
        else
        {
          if( $head_splicer eq 'AG' )
          {
            $head=1;
          }
          else
          {
            $head=-1;
          }
        }

        #Check tail splicer
        if( $exon_number == $exon_count_max )
        {
          $tail=0;
        }
        else
        {
          if( $tail_splicer eq 'GT' )
          {
            $tail=1;
          }  
          else
          {
            $tail=-1;
          }
        }

        if ($sth3 = $dbh->prepare("insert tb_splicer (head,tail,exon_id,datasource_id) value ($head,$tail,$exon_id,$datasource_id)"))
        {
          if ($sth3->execute())
          {
            $db_count++;
            print "Writen db records: $db_count\r";
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
      }
    }
    else
    {
      print "\nError in:$_\n";
      print "Can't execute SQL statement: $DBI::errstr\n";
      next;
    }
  }
  else
  {
    print "\nCan't prepare SQL statement: $DBI::errstr\n";
    next;
  }
}
print "\n";

$sth3->finish;
$sth2->finish;
$sth1->finish;

$dbh->disconnect;

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
}

