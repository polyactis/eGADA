#Make a base image without eGADA

all:
	date
	ls bashrc
	docker login
	docker pull ubuntu:20.04
	docker build -t polyactis/egada:u20base .
	docker push polyactis/egada:u20base
	date
