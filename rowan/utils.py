import os

import rowan


def get_api_key() -> str:
    api_key = os.environ.get("ROWAN_API_KEY")
    if api_key is not None:
        return api_key
    elif hasattr(rowan, "api_key"):
        return rowan.api_key
    else:
        raise ValueError(
            "No API key provided. You can set your API key using 'rowan.api_key = <API-KEY>', or you can set the environment variable ROWAN_API_KEY=<API-KEY>)."
        )
