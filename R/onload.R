'.onAttach' <- function(lib, pkg="goftte")
  {    
    desc <- packageDescription(pkg)
    packageStartupMessage("Loading '", desc$Package, "' version ",desc$Version);
  }
