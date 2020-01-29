test:
	nextflow run main.nf -profile docker,test

docker_push:
	sudo docker login
	docker push nfcore/kmermaid:dev

docker_build:
	@docker build \
		--build-arg VCS_REF=`git rev-parse --short HEAD`  \
		--build-arg BUILD_DATE=`date -u +"%Y-%m-%dT%H:%M:%SZ"` \
		-t nfcore/kmermaid:dev .
