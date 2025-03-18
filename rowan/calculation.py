from typing import Any

import stjames

from .utils import api_client


class Calculation:
    @classmethod
    def retrieve(cls, uuid: stjames.UUID) -> dict[str, Any]:
        with api_client() as client:
            response = client.get(f"/calculation/{uuid}/stjames")
            response.raise_for_status()
            return response.json()
