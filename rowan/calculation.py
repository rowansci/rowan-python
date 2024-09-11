import stjames

from .utils import api_client

class Calculation:
    @classmethod
    def retrieve(cls, uuid: stjames.UUID) -> dict:
        with api_client() as client:
            response = client.get(f"/calculation/{uuid}/stjames")
            response.raise_for_status()
            return response.json()
