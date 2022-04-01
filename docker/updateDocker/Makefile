# Add Accucopy to the base image

tag_name = v1.0.1

all:
	date
	docker login
	docker login docker.yfish.x
	-mkdir binary
	cd ../../src/; make;
	rm -rf ./binary/*
	cd ../../src/ && cp -apr eGADA eGADA.so testGADA.py ../docker/updateDocker/binary/
	docker pull polyactis/egada:u20base
	docker build -t polyactis/egada:latest .
	docker tag polyactis/egada:latest polyactis/egada:$(tag_name)
	docker push polyactis/egada:$(tag_name)
	docker push polyactis/egada:latest
	docker tag polyactis/egada:latest docker.yfish.x/egada:$(tag_name)
	docker tag docker.yfish.x/egada:$(tag_name) docker.yfish.x/egada:latest
	docker push docker.yfish.x/egada
	date
