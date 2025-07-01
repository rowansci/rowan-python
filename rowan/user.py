from typing import Any

from .utils import api_client


class User:
    """
    User class
    """

    def __init__(
        self,
        uuid: str,
        username: str,
        email: str,
        firstname: str | None = None,
        lastname: str | None = None,
        weekly_credits: float | None = None,
        credits: float | None = None,
        billing_name: str | None = None,
        billing_address: str | None = None,
        credit_balance_warning: float | None = None,
        organization: dict[str, Any] | None = None,
        organization_role: dict[str, Any] | None = None,
        individual_subscription: dict[str, Any] | None = None,
        **kwargs,
    ):
        """
        Initialize a Folder object.

        :param uuid: Unique identifier for the folder.
        :type uuid: UUID
        :param name: Name of the folder.
        :type name: str or None
        :param parent_uuid: UUID of the parent folder.
        :type parent_uuid: UUID or None
        :param notes: Optional notes about the folder.
        :type notes: str
        :param starred: Whether the folder is starred.
        :type starred: bool
        :param public: Whether the folder is public.
        :type public: bool
        """
        self.uuid = uuid
        self.username = username
        self.name = f"{firstname} {lastname}"
        self.email = email
        self.weekly_credits = weekly_credits
        self.credits = credits
        self.billing_name = billing_name
        self.billing_address = billing_address
        self.credit_balance_warning = credit_balance_warning
        self.organization = organization["name"] if organization is not None else None
        self.organization_weekly_credits = (
            organization["weekly_credits"] if organization is not None else None
        )
        self.organization_credits = organization["credits"] if organization is not None else None
        self.organization_role = (
            organization_role["name"] if organization_role is not None else None
        )
        self.individual_subscription = (
            individual_subscription["subscription_plan"]["name"]
            if individual_subscription is not None
            else None
        )

    def __repr__(self):
        return f"User(uuid={self.uuid}, username={self.username}, email={self.email})"

    def credits_available_string(self):
        individual_credits = f"Weekly Credits: {self.weekly_credits}\nCredits: {self.credits}"
        if self.organization is not None:
            print(
                individual_credits
                + f"\nOrganization Weekly Credits: {self.organization_weekly_credits}\n"
                f"Organization Credits: {self.organization_credits}"
            )
            return
        print(individual_credits)


def me():
    """
    Returns the current user
    """

    with api_client() as client:
        response = client.get("/user/me")
        response.raise_for_status()

    return User(**response.json())
