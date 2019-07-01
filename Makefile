test:
	nextflow run main.nf -profile docker,test

docker_push:
	sudo docker login
	docker push czbiohub/nf-kmer-similarity:dev

docker_build:
	@docker build \
		--build-arg VCS_REF=`git rev-parse --short HEAD`  \
		--build-arg BUILD_DATE=`date -u +"%Y-%m-%dT%H:%M:%SZ"` \
<<<<<<< HEAD
		-t czbiohub/nf-kmer-similarity:dev .
=======
		-t czbiohub/nf-kmer-similarity .

docker_push:
	sudo docker login
	sudo docker push czbiohub/nf-kmer-similarity:olgabot-dayhoff
	docker images
>>>>>>> Update dockerfile to use olgabot/dayhoff branch"
