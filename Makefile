RCMD := Rscript -e
.PHONY: install
install: ## Install all HiC* packages and dependencies with pak.
	@echo "🚀 Installing HiCExperiment package"
	$(RCMD) 'install.packages("pak", repos = "https://r-lib.github.io/p/pak/devel/")'
	$(RCMD) 'pak::pkg_install(".", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/HiCool", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/HiContacts", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/HiContactsData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/fourDNData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/DNAZooData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'

.PHONY: check
check: deps ## Check if HiCExperiment package can be installed.
	@echo "🔍 Checking HiCExperiment"
	$(RCMD) 'devtools::check()'

deps: ## Install all HiCExperiment dependencies
	@echo "🔗 Installing dependencies"
	$(RCMD) 'devtools::install_dev_deps(".")'

.PHONY: docker-build
docker-build: ## Build the Docker image.
	@echo "🐋 Building docker image"
	docker build -t r-hicexperiment .

.PHONY: test
test: install ## Test the code with testthat.
	@echo "🧪 Testing code: Running testthat"
	$(RCMD) 'devtools::test()'

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
