#include "cartmesh/cutcell/CutCellMeshFingerprint.hpp"

#include <bit>
#include <iomanip>
#include <sstream>

namespace cartmesh {
namespace {

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

void hash_u64(std::uint64_t& hash, std::uint64_t value) noexcept {
    for (std::uint32_t shift = 0; shift < 64U; shift += 8U) {
        hash_byte(hash, static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
}

void hash_double(std::uint64_t& hash, double value) noexcept {
    hash_u64(hash, std::bit_cast<std::uint64_t>(value));
}

void hash_vec3(std::uint64_t& hash, const Vec3& value) noexcept {
    hash_double(hash, value.x);
    hash_double(hash, value.y);
    hash_double(hash, value.z);
}

void hash_mesh_payload(std::uint64_t& hash,
                       const ConvexCutCellMesh& mesh) noexcept {
    for (const auto& cell : mesh.fluid_cells) {
        hash_u64(hash, cell.background_cell_id);
        hash_double(hash, cell.volume);
        hash_double(hash, cell.volume_fraction);
        hash_vec3(hash, cell.centroid);
        hash_u64(hash, cell.fluid_component_count);
        for (const auto& face : cell.cartesian_faces) {
            hash_double(hash, face.area);
            hash_vec3(hash, face.centroid);
            hash_vec3(hash, face.outward_normal);
            for (const auto& loop : face.oriented_boundary_loops) {
                hash_u64(hash, loop.size());
                for (const auto& vertex : loop) hash_vec3(hash, vertex);
            }
        }
        for (const auto& face : cell.embedded_boundary_faces) {
            hash_u64(hash, face.boundary_id);
            hash_double(hash, face.area);
            hash_vec3(hash, face.centroid);
            hash_vec3(hash, face.outward_normal);
            for (const auto& vertex : face.vertices) hash_vec3(hash, vertex);
        }
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            hash_u64(hash, piece.component_id);
            hash_u64(hash, piece.global_region_id);
            hash_double(hash, piece.geometry.volume);
            for (const auto& vertex : piece.polyhedron.vertices) {
                hash_vec3(hash, vertex);
            }
            for (const auto& face : piece.polyhedron.faces) {
                hash_byte(hash, static_cast<std::uint8_t>(face.kind));
                hash_u64(hash, face.source_id);
                for (const auto vertex : face.vertex_indices) {
                    hash_u64(hash, vertex);
                }
            }
        }
        for (const auto region_id : cell.fluid_component_region_ids) {
            hash_u64(hash, region_id);
        }
    }
    for (const auto& face : mesh.internal_faces) {
        hash_u64(hash, face.first_background_cell_id);
        hash_u64(hash, face.second_background_cell_id);
        hash_double(hash, face.area);
        hash_vec3(hash, face.centroid);
        hash_vec3(hash, face.normal);
    }
    for (const auto& face : mesh.component_internal_faces) {
        hash_u64(hash, face.first_background_cell_id);
        hash_u64(hash, face.second_background_cell_id);
        hash_u64(hash, face.first_component_id);
        hash_u64(hash, face.second_component_id);
        hash_u64(hash, face.global_region_id);
        hash_double(hash, face.area);
        hash_vec3(hash, face.centroid);
    }
    for (const auto volume : mesh.global_fluid_region_volumes) {
        hash_double(hash, volume);
    }
}

} // namespace

std::uint64_t cut_cell_mesh_fingerprint_fnv1a64(
    const LinearOctree& tree, const ConvexCutCellMesh& mesh) noexcept {
    std::uint64_t hash = 14695981039346656037ULL;
    hash_vec3(hash, tree.domain().minimum());
    hash_vec3(hash, tree.domain().maximum());
    hash_u64(hash, tree.base_level());
    hash_u64(hash, tree.maximum_level());
    for (const auto code : tree.leaf_codes()) hash_u64(hash, code);
    hash_mesh_payload(hash, mesh);
    return hash;
}

std::string fingerprint_hex(std::uint64_t fingerprint) {
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(16) << fingerprint;
    return output.str();
}

} // namespace cartmesh
