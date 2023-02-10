RCMD := Rscript -e
.PHONY: install
install: deps ## Install package and dependencies.
	@echo "🚀 Installing package"
	$(RCMD) 'devtools::install()'

.PHONY: check
check: deps ## Check if package can be installed.
	@echo "🔍 Checking HiCExperiment"
	$(RCMD) 'devtools::check()'

deps: ## Install all dependencies
	@echo "🔗 Installing dependencies"
	$(RCMD) 'devtools::install_dev_deps(".")'

.PHONY: docker-build
docker-build: ## Build the Docker image.
	@echo "🐋 Building docker image"
	docker build -t r-hicexperiment .

.PHONY: test
test: install ## Test the code with testthat.
	@echo "🧪 Testing code: Running testthat"
	$(RCMD) 'devtools::load_all();devtools::test()'

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
