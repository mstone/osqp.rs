{
  description = "c2rust port of the Operator Splitting Quadratic Programming (OSQP) solver";

  inputs.c2rust.url = "github:mstone/c2rust";

  inputs.crane.url = "github:ipetkov/crane";
  inputs.crane.inputs.flake-utils.follows = "flake-utils";
  inputs.crane.inputs.nixpkgs.follows = "nixpkgs";

  inputs.flake-utils.url = "github:numtide/flake-utils";

  inputs.nixpkgs.url = "nixpkgs/nixpkgs-unstable";

  inputs.rust-overlay.url = "github:oxalica/rust-overlay";
  inputs.rust-overlay.inputs.flake-utils.follows = "flake-utils";
  inputs.rust-overlay.inputs.nixpkgs.follows = "nixpkgs";

  outputs = {self, nixpkgs, c2rust, crane, rust-overlay, flake-utils}:
    flake-utils.lib.simpleFlake {
      inherit self nixpkgs;
      name = "osqp";
      systems = flake-utils.lib.allSystems;
      preOverlays = [ rust-overlay.overlay ];
      overlay = final: prev: {
        osqp = rec {
          osqp = lib.osqp { isShell = false; };
          devShell = lib.osqp { isShell = true; };
          defaultPackage = osqp;

          rust = with final; with pkgs; (rust-bin.nightly.latest.minimal);

          lib.osqp = { isShell, subpkg ? "osqp", subdir ? "." }:
            let
              buildInputs = with final; with pkgs; [
                rust
                cmake
                c2rust.legacyPackages."${final.system}".defaultPackage
              ] ++ final.lib.optionals isShell [
                entr
              ] ++ final.lib.optionals stdenv.isDarwin (with darwin.apple_sdk.frameworks; [
              ]) ++ final.lib.optionals stdenv.isLinux ([
              ]);
            in with final; with pkgs; crane.lib.${final.system}.buildPackage {
            pname = "${subpkg}";
            version = "0.6.2";

            src = self;

            inherit buildInputs;
            dontUseCmakeConfigure = true;
            doCheck = false;
          };
        };
      };
    };
}
