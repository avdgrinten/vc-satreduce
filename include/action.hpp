#pragma once

namespace vc_bnb {

enum action_type {
    invalid,
    add_edge,
    pick,
    unpick,
    unify,
    unified,
    funnel,
    constraint_created,
    constraint_removed,
    removed_from_constraint_pick,
    removed_from_constraint_unpick,
    constraint_disabled,
};

} // namespace vc_bnb
