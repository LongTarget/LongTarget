#!/usr/bin/perl

use strict;
use Getopt::Std;
use Time::Local;
use DBI;
use Bio::EnsEMBL::Registry;

getopts('h:u:p:d:f:s:c:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $db_file, $datasource_id, $species_id);

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

if ( defined($opts{f}) )
{
  $db_file=$opts{f};
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
  $datasource_id=2;#lncrna db
}

if ( defined($opts{c}))
{
  $species_id=$opts{c};
}
else
{
  $species_id=1;
}

if (!open(DBFILE, $db_file))
{
  die "Can't open $db_file.\n";
}

my $dbh = DBI->connect( "dbi:mysql:dbname=$dbname;host=$dbhost", $dbuser, $dbpasswd);
my $sth1, $sth2;
if ( !defined $dbh )
{
  close DBFILE;
  die "Cannot connect to database  $dbhost $dbname $dbuser.\n";
}

my $db_count=0;
my $last_gene_name='';
my $tb_gene_name='';
my $gene_id=0;
my $do_alias=0;
my @fields;
while (<DBFILE>)
{
  if ( /^#/ || /^\s+/ )
  {
    next;
  }

  chomp;
  my ($gene_name, $alias, $attr1, $attr2) = split /\t/;

  my $sqlstr;

  if ($gene_name ne $last_gene_name)
  {
    $last_gene_name=$gene_name;
    @fields=split /,\s*/, $alias;
    push(@fields, $gene_name);
    my $rec;
    foreach $rec (@fields)
    {
      $sqlstr="select gene_id from tb_gene where gene_name=\"$rec\"";
      $sth1=$dbh->prepare($sqlstr)
        or die "Can't prepare SQL statement: $DBI::errstr\n$sqlstr\n";
      $sth1->execute();
      if(($gene_id) = $sth1->fetchrow_array())
      {
      	$tb_gene_name=$rec;
      	$do_alias=1;
      	last;
      }
    }
    if( !defined($gene_id))
    {
      print "Can't find $gene_name, $alias\n";
    }
    else
    {
      foreach $rec (@fields)
      {
      	if( $rec eq $tb_gene_name)
      	{
      	  next;
      	}
        if ($sth2 = $dbh->prepare("insert tb_alias (gene_id,alias_name,datasource_id) value ($gene_id,\"$rec\",$datasource_id)"))
        {
          if ($sth2->execute())
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
}
$sth2->finish;
$sth1->finish;

$dbh->disconnect;

close(DBFILE);

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
  print "-f [lncrnadb file name], default is '0'.\n";
  print "-o [Output result path].\n";
}

