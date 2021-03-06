#define CATCH_CONFIG_MAIN

#include <catch.hpp>
#include <GkylRange.h>
#include <cmath>
#include <vector>

TEST_CASE("Tests for 1D range object", "[range1D]") {
  GkylRange_t crange;
  crange.ndim = 1;
  crange.lower[0] = 1; crange.upper[0] = 5;
  Gkyl::calcIndexerCoeff(crange);

  Gkyl::Range range(&crange);

  REQUIRE( range.ndim() == 1 );
  for (unsigned i=0; i<range.ndim(); ++i) {
    REQUIRE( range.lower(i) == crange.lower[i] );
    REQUIRE( range.upper(i) == crange.upper[i] );
    REQUIRE( range.shape(i) == (crange.upper[i]-crange.lower[i]+1) );
  }
  REQUIRE( range.volume() == 5 );

  Gkyl::Indexer<1> idxr(&crange);
  Gkyl::GenIndexer genIdxr(&crange);

  int idx[1];
  int count = 0;
  for (int i=range.lower(0); i<=range.upper(0); ++i) {
    REQUIRE( idxr.index(i) == count );

    idx[0] = i;
    REQUIRE( genIdxr.index(idx) == count );
    
    count++;
  }

  for (int i=range.lower(0); i<=range.upper(0); ++i) {
    idxr.invIndex(idxr.index(i), idx);
    REQUIRE( idx[0] == i );
  }
      
}

TEST_CASE("Tests for 2D range object row-major", "[range2D-row-major]") {
  GkylRange_t crange;
  crange.ndim = 2;
  crange.lower[0] = 1; crange.upper[0] = 5;
  crange.lower[1] = 2; crange.upper[1] = 10;
  Gkyl::calcIndexerCoeff(crange);

  Gkyl::Range range(&crange);
  Gkyl::GenIndexer genIdxr(&crange);  

  REQUIRE( range.ndim() == 2 );
  for (unsigned i=0; i<range.ndim(); ++i) {
    REQUIRE( range.lower(i) == crange.lower[i] );
    REQUIRE( range.upper(i) == crange.upper[i] );
    REQUIRE( range.shape(i) == (crange.upper[i]-crange.lower[i]+1) );
  }
  REQUIRE( range.volume() == 5*9 );

  Gkyl::Indexer<2> idxr(&crange);

  int idx[2];
  int count = 0;
  for (int i=range.lower(0); i<=range.upper(0); ++i)
    for (int j=range.lower(1); j<=range.upper(1); ++j) {
      
      REQUIRE( idxr.index(i,j) == count );

      idx[0] = i; idx[1] = j;
      REQUIRE( genIdxr.index(idx) == count );
      
      count++;
    }

  for (int i=range.lower(0); i<=range.upper(0); ++i)
    for (int j=range.lower(1); j<=range.upper(1); ++j) {

      idxr.invIndex(idxr.index(i,j), idx);

      REQUIRE( idx[0] == i );
      REQUIRE( idx[1] == j );

      idx[0] = 0; idx[1] = 0;
      genIdxr.invIndex(idxr.index(i,j), idx);

      REQUIRE( idx[0] == i );
      REQUIRE( idx[1] == j );
    }
}

TEST_CASE("Tests for 2D range object col-major", "[range2D-col-major]") {
  GkylRange_t crange;
  crange.ndim = 2;
  crange.lower[0] = 1; crange.upper[0] = 5;
  crange.lower[1] = 2; crange.upper[1] = 10;
  Gkyl::calcIndexerCoeff(crange);

  Gkyl::Range range(&crange);
  Gkyl::GenIndexer genIdxr(&crange, Gkyl::Layout::colMajor);

  REQUIRE( range.ndim() == 2 );
  for (unsigned i=0; i<range.ndim(); ++i) {
    REQUIRE( range.lower(i) == crange.lower[i] );
    REQUIRE( range.upper(i) == crange.upper[i] );
    REQUIRE( range.shape(i) == (crange.upper[i]-crange.lower[i]+1) );
  }
  REQUIRE( range.volume() == 5*9 );

  Gkyl::Indexer<2> idxr(&crange, Gkyl::Layout::colMajor);

  int idx[2];
  int count = 0;
  for (int j=range.lower(1); j<=range.upper(1); ++j)  
    for (int i=range.lower(0); i<=range.upper(0); ++i) {
      REQUIRE( idxr.index(i,j) == count );

      idx[0] = i; idx[1] = j;
      REQUIRE( genIdxr.index(idx) == count );
      
      count++;
    }

  for (int j=range.lower(1); j<=range.upper(1); ++j)  
    for (int i=range.lower(0); i<=range.upper(0); ++i) {

      idxr.invIndex(idxr.index(i,j), idx);
      
      REQUIRE( idx[0] == i );
      REQUIRE( idx[1] == j );
    }
      
}

TEST_CASE("Tests for 3D range object row-major", "[range3D-row-major]") {
  GkylRange_t crange;
  crange.ndim = 3;
  crange.lower[0] = 1; crange.upper[0] = 5;
  crange.lower[1] = 2; crange.upper[1] = 10;
  crange.lower[2] = -1; crange.upper[2] = 4;
  Gkyl::calcIndexerCoeff(crange);

  Gkyl::Range range(&crange);
  Gkyl::GenIndexer genIdxr(&crange);  

  REQUIRE( range.ndim() == 3 );
  for (unsigned i=0; i<range.ndim(); ++i) {
    REQUIRE( range.lower(i) == crange.lower[i] );
    REQUIRE( range.upper(i) == crange.upper[i] );
    REQUIRE( range.shape(i) == (crange.upper[i]-crange.lower[i]+1) );
  }
  REQUIRE( range.volume() == 5*9*6 );

  Gkyl::Indexer<3> idxr(&crange);

  int idx[3];
  int count = 0;
  for (int i=range.lower(0); i<=range.upper(0); ++i)
    for (int j=range.lower(1); j<=range.upper(1); ++j)
      for (int k=range.lower(2); k<=range.upper(2); ++k) {
      
        REQUIRE( idxr.index(i,j,k) == count );

        idx[0] = i; idx[1] = j; idx[2] = k;
        REQUIRE( genIdxr.index(idx) == count );
      
        count++;
      }

  for (int i=range.lower(0); i<=range.upper(0); ++i)
    for (int j=range.lower(1); j<=range.upper(1); ++j)
      for (int k=range.lower(2); k<=range.upper(2); ++k) {

        idxr.invIndex(idxr.index(i,j,k), idx);
        
        REQUIRE( idx[0] == i );
        REQUIRE( idx[1] == j );
        REQUIRE( idx[2] == k );

        idx[0] = 0; idx[1] = 0; idx[2] = 0;
        genIdxr.invIndex(idxr.index(i,j,k), idx);
        
        REQUIRE( idx[0] == i );
        REQUIRE( idx[1] == j );
        REQUIRE( idx[2] == k );
      }
      
}

TEST_CASE("Tests for 3D range object col-major", "[range3D-col-major]") {
  GkylRange_t crange;
  crange.ndim = 3;
  crange.lower[0] = 1; crange.upper[0] = 5;
  crange.lower[1] = 2; crange.upper[1] = 10;
  crange.lower[2] = -1; crange.upper[2] = 4;
  Gkyl::calcIndexerCoeff(crange);

  Gkyl::Range range(&crange);
  Gkyl::GenIndexer genIdxr = range.genIndexer(Gkyl::Layout::colMajor);

  REQUIRE( range.ndim() == 3 );
  for (unsigned i=0; i<range.ndim(); ++i) {
    REQUIRE( range.lower(i) == crange.lower[i] );
    REQUIRE( range.upper(i) == crange.upper[i] );
    REQUIRE( range.shape(i) == (crange.upper[i]-crange.lower[i]+1) );
  }
  REQUIRE( range.volume() == 5*9*6 );

  Gkyl::Indexer<3> idxr(&crange, Gkyl::Layout::colMajor);

  int idx[3];
  int count = 0;

  for (int k=range.lower(2); k<=range.upper(2); ++k)    
    for (int j=range.lower(1); j<=range.upper(1); ++j)
      for (int i=range.lower(0); i<=range.upper(0); ++i) {
      
        REQUIRE( idxr.index(i,j,k) == count );

        idx[0] = i; idx[1] = j; idx[2] = k;
        REQUIRE( genIdxr.index(idx) == count );
      
        count++;
      }

  for (int i=range.lower(0); i<=range.upper(0); ++i)
    for (int j=range.lower(1); j<=range.upper(1); ++j)
      for (int k=range.lower(2); k<=range.upper(2); ++k) {

        idxr.invIndex(idxr.index(i,j,k), idx);
        
        REQUIRE( idx[0] == i );
        REQUIRE( idx[1] == j );
        REQUIRE( idx[2] == k );

        idx[0] = 0; idx[1] = 0; idx[2] = 0;
        genIdxr.invIndex(idxr.index(i,j,k), idx);
        
        REQUIRE( idx[0] == i );
        REQUIRE( idx[1] == j );
        REQUIRE( idx[2] == k );
      }
      
}

