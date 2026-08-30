#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"
#include "cartmesh2d/hybrid/TransitionCanonicalization2D.hpp"
#include "cartmesh2d/quality/QualityContract2D.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cartmesh2d;
namespace {
int failures=0;
void check(bool condition,const char* message) {
    if (!condition) {++failures;std::cerr<<"FAIL: "<<message<<'\n';}
}
bool equal(Point2D a,Point2D b) {return a.x==b.x && a.y==b.y;}
}

int main() {
    check(TransitionCanonicalizationPolicy2D{}.minimumFaceFraction==
          QualityContract2D{}.transition.faceOverLocalBackgroundH.hard,
          "sampling criterion equals the unchanged Q1 face/local_h hard limit");
    {
        const FeatureOwner2D ownerA{ConstructionSourceKind2D::WallSegment, 7U, 1U};
        const FeatureOwner2D ownerB{ConstructionSourceKind2D::WallSegment, 8U, 1U};
        check(constructionFeaturesCompatible(
                  ConstructionFeatureClass2D::ConvexSharp, ownerA,
                  ConstructionFeatureClass2D::ConvexSharp, ownerA).compatible,
              "same sharp feature owner is compatible");
        check(!constructionFeaturesCompatible(
                   ConstructionFeatureClass2D::ConvexSharp, ownerA,
                   ConstructionFeatureClass2D::ConcaveSharp, ownerA).compatible &&
              !constructionFeaturesCompatible(
                   ConstructionFeatureClass2D::ConcaveSharp, ownerA,
                   ConstructionFeatureClass2D::ConcaveSharp, ownerB).compatible,
              "sharp/concave type and owner conflicts are explicit");
        check(!constructionFeaturesCompatible(
                   ConstructionFeatureClass2D::GapSideA, ownerA,
                   ConstructionFeatureClass2D::GapSideB, ownerA).compatible &&
              !constructionFeaturesCompatible(
                   ConstructionFeatureClass2D::GapSideA, ownerA,
                   ConstructionFeatureClass2D::GapSideA, ownerB).compatible &&
              constructionFeaturesCompatible(
                   ConstructionFeatureClass2D::GapSideA, ownerA,
                   ConstructionFeatureClass2D::GapSideA, ownerA).compatible,
              "gap sides require identical side class and owner");
        check(constructionFeaturesCompatible(
                  ConstructionFeatureClass2D::Smooth, std::nullopt,
                  ConstructionFeatureClass2D::ConvexSharp, ownerA).compatible,
              "a movable sample may reuse an immutable sharp anchor");
    }
    {
        ConstructionVertexStore2D store;
        const auto left = store.addShadowVertex(
            {-0.25, 0.0}, 1.0, ConstructionFeatureClass2D::Smooth,
            std::nullopt, 3U);
        const auto right = store.addShadowVertex(
            {0.25, 0.0}, 1.0, ConstructionFeatureClass2D::Smooth,
            std::nullopt, 3U);
        const StableVertexKey2D key{StableVertexKeyKind2D::GridVertex,
                                    0U, 0U, 4U, 5U};
        store.bindExactKey(key, left);
        store.bindExactKey(key, left);
        check(store.resolveExactKey(key) == left,
              "typed exact key binding is idempotent");
        bool conflict = false;
        try { store.bindExactKey(key, right); }
        catch (const std::runtime_error&) { conflict = true; }
        check(conflict, "typed exact key conflict fails closed");
        const auto decision = store.decideProximity(
            {{0.0, 0.0}, 1.0, 0.3, ConstructionFeatureClass2D::Smooth,
             std::nullopt, 3U, false});
        check(decision.canonicalId == left && equal(decision.canonicalPoint, {-0.25, 0.0}),
              "equidistant proximity decision is geometric and deterministic");
    }
    {
        ConstructionVertexStore2D store;
        const FeatureOwner2D sideA{ConstructionSourceKind2D::WallSegment, 10U, 0U};
        const FeatureOwner2D sideB{ConstructionSourceKind2D::WallSegment, 11U, 0U};
        (void)store.addShadowVertex({0.0, 0.0}, 1.0,
                                    ConstructionFeatureClass2D::GapSideA,
                                    sideA, 9U);
        const auto gapDecision = store.decideProximity(
            {{1.e-3, 0.0}, 1.0, 0.1, ConstructionFeatureClass2D::GapSideB,
             sideB, 9U, false});
        check(!gapDecision.canonicalId &&
              equal(gapDecision.canonicalPoint, {1.e-3, 0.0}),
              "opposite gap sides remain distinct inside snap radius");
        const auto sharpDecision = store.decideProximity(
            {{1.e-3, 0.0}, 1.0, 0.1, ConstructionFeatureClass2D::ConvexSharp,
             sideA, 9U, true});
        check(!sharpDecision.canonicalId,
              "immutable sharp proposal is never displaced");
    }
    for (const double scale:{1.e-6,1.0,1.e6}) {
        IntersectionRegistry2D registry({1.e-6});
        // Exact Q1 failure coordinates, not a zeroed/truncated fixture.
        const Point2D intersection{-2.326100423965795*scale,9.769353911109006e-9*scale};
        const Point2D anchor{-2.3261004232233007*scale,0.0};
        const double h=.0875*scale;
        (void)registry.addCanonicalVertex(anchor,h,IntersectionFeature2D::TransitionVertex,
                                         std::nullopt,7U);
        check(equal(registry.canonicalize(intersection,h,
              IntersectionSource2D::TransitionEnvelopeCartesian,363U,
              IntersectionFeature2D::Smooth,std::nullopt,7U),anchor),
              "the micro-edge endpoint resolves at every physical scale");
        check(registry.records().size()==1U,"one changed intersection is retained");
        const auto& record=registry.records().front();
        check(equal(record.originalPoint,intersection) && record.localH==h &&
              record.sourceId==363U && record.supportId==7U &&
              record.feature==IntersectionFeature2D::Smooth,
              "complete original source, position, local_h and feature provenance");
        check(std::abs(record.displacement/h-9.79752897103936e-9/.0875)<1.e-14,
              "displacement/local_h is scale invariant");
        check(equal(registry.canonicalize(intersection,h,
              IntersectionSource2D::WallCartesian,1U,IntersectionFeature2D::Smooth,
              std::nullopt,8U),intersection),"nearby nonincident support cannot snap");
        check(equal(registry.canonicalize(intersection,h,
              IntersectionSource2D::WallCartesian,1U,IntersectionFeature2D::WallSharpCorner,
              3U,7U),intersection),"a real sharp corner cannot be moved");
        check(equal(registry.canonicalize(intersection,h,
              IntersectionSource2D::WallCartesian,1U,IntersectionFeature2D::WallConcaveCorner,
              4U,7U),intersection),"a concave feature cannot be moved");
        const auto repeat=intersectionRecordsToJson(registry.records());
        check(repeat==intersectionRecordsToJson(registry.records()),"provenance serialization is deterministic");
        check(repeat.find("displacement_over_local_h")!=std::string::npos,
              "machine-readable provenance includes dimensionless displacement");
        const auto& shadow=registry.shadowVertexStore();
        check(shadow.records().size()==registry.vertices().size(),
              "R1-A shadow store has one stable record per canonical vertex");
        check(shadow.indexProfile().queryCount>0U &&
              shadow.indexProfile().maximumQueryCandidateCount<=registry.vertices().size(),
              "R1-A feature index shadows every legacy canonicalization query");
        for (std::size_t i=0;i<shadow.records().size();++i) {
            check(shadow.records()[i].id==i &&
                  shadow.records()[i].key.kind==StableVertexKeyKind2D::LegacyCanonical,
                  "shadow stable ids and typed legacy keys are deterministic");
        }
    }
    {
        IntersectionRegistry2D registry({.01});
        (void)registry.addCanonicalVertex({0,0},1.e-6,IntersectionFeature2D::CartesianGridVertex);
        check(equal(registry.canonicalize({1.e-5,0},1.,IntersectionSource2D::WallCartesian,0),{1.e-5,0}),
              "a coarse side cannot absorb a fine resolved feature");
        bool threw=false;
        try {(void)registry.canonicalize({0,0},0.,IntersectionSource2D::WallCartesian,0);}
        catch (const std::invalid_argument&) {threw=true;}
        check(threw,"zero local_h fails explicitly");
    }
    {
        const Point2D near{9.e-9,.8};
        const Domain2D domain{{{-1,-1},{1,1}}};
        const BoundaryRegion2D fixed(BoundaryLoop({{-.2,-.2},{.2,-.2},{.2,.2},{-.2,.2}}));
        std::vector<BoundaryLoop> front{BoundaryLoop({{-.8,-.8},{.8,-.8},{.8,.8},near,{-.8,.8}})};
        std::vector<Polygon2D> cells{
            {{{-.8,-.8},{.8,-.8},{.2,-.2},{-.2,-.2}}},
            {{{.8,-.8},{.8,.8},{.2,.2},{.2,-.2}}},
            {{{.8,.8},near,{-.8,.8},{-.2,.2},{.2,.2}}},
            {{{-.8,.8},{-.8,-.8},{-.2,-.2},{-.2,.2}}}};
        double before=0;for (const auto& cell:cells) before+=cell.area();
        IntersectionRegistry2D registry({.01});
        std::string error;
        check(canonicalizeTransitionEnvelope2D(front,cells,fixed,domain,2U,.5,false,
              registry,{},error),"minimal smooth envelope can be resampled before cutting");
        check(equal(front.front().vertices()[3],{0,.8}),"near-axis transition sample uses canonical grid line");
        check(equal(cells[2].vertices[1],{0,.8}),"the adjacent transition polygon uses the identical canonical point");
        double after=0;for (const auto& cell:cells) after+=cell.area();
        check(std::abs(after-before)<1.e-14 && cells.size()==4U,"common sampling conserves area and does not delete cells");
        check(equal(front.front().vertices()[2],{.8,.8}),"sharp envelope corner is preserved");
    }
    return failures==0?0:1;
}
