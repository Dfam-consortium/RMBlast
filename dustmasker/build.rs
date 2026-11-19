// Embeds the git commit in --version; see build-support/git_build_id.rs.
include!("../build-support/git_build_id.rs");

fn main() {
    emit_build_id();
}
