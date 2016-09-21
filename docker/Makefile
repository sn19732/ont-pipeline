IMAGE=ont-assembly-polish

build: Dockerfile
	docker build --pull -t $(IMAGE) .

.PHONY: push
push:
	docker tag $(IMAGE):latest docker-registry.oxfordnanolabs.local:5000/$(IMAGE):latest
	docker push docker-registry.oxfordnanolabs.local:5000/$(IMAGE)
